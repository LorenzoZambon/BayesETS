##############################################################################
### Adaptive Multiple Importance Sampling with MVT Mixture ###

# Initialize a mixture of K MVT proposals from a base single-component proposal
# base_params: list(mus, Sigma, df) from init_joint_params
# Returns: list(K, weights, mus_list, Sigma_list, df)
init_mixture_proposal <- function(base_params, K = 3, jitter_scale = 0.5) {
  d <- length(base_params$mus)
  weights <- rep(1 / K, K)
  mus_list <- vector("list", K)
  Sigma_list <- vector("list", K)
  param_names <- names(base_params$mus)

  mus_list[[1]] <- base_params$mus
  Sigma_list[[1]] <- base_params$Sigma

  if (K > 1) {
    sd_vec <- sqrt(pmax(diag(base_params$Sigma), 1e-6))
    for (k in 2:K) {
      mus_list[[k]] <- base_params$mus + stats::rnorm(d) * sd_vec * jitter_scale
      Sigma_list[[k]] <- base_params$Sigma * (1 + 0.5 * (k - 1) / K)
    }
  }

  for (k in seq_len(K)) {
    names(mus_list[[k]]) <- param_names
    rownames(Sigma_list[[k]]) <- colnames(Sigma_list[[k]]) <- param_names
  }

  list(
    K = K, weights = weights,
    mus_list = mus_list, Sigma_list = Sigma_list,
    df = base_params$df
  )
}


# Evaluate log-density of a MVT mixture in unconstrained space (vectorized)
# x: matrix (N x d), mixture_params: mixture specification
# Returns: numeric vector of length N
log_density_mixture_unc <- function(x, mixture_params) {
  K <- mixture_params$K
  df <- mixture_params$df

  log_comp <- matrix(NA_real_, nrow = nrow(x), ncol = K)
  for (k in seq_len(K)) {
    log_comp[, k] <- log(mixture_params$weights[k]) +
      mvtnorm::dmvt(
        x,
        delta = mixture_params$mus_list[[k]],
        sigma = mixture_params$Sigma_list[[k]],
        df = df, log = TRUE
      )
  }

  # Log-sum-exp across components
  max_log <- apply(log_comp, 1, max)
  max_log + log(rowSums(exp(log_comp - max_log)))
}


# Sample N particles from the MVT mixture and transform to constrained space
# Returns same structure as draw_from_joint_proposal plus log_density_unc and log_jac
draw_from_mixture_proposal <- function(N, mixture_params, theta_names,
                                       phi_min, phi_max) {
  K <- mixture_params$K
  df <- mixture_params$df
  param_names <- names(mixture_params$mus_list[[1]])
  d <- length(param_names)

  # Assign each particle to a component
  comp_assignment <- sample(K, size = N, replace = TRUE,
                            prob = mixture_params$weights)
  n_per_comp <- tabulate(comp_assignment, nbins = K)

  samps_unc <- matrix(NA_real_, nrow = N, ncol = d)
  colnames(samps_unc) <- param_names

  for (k in seq_len(K)) {
    if (n_per_comp[k] > 0) {
      idx <- which(comp_assignment == k)
      samps_unc[idx, ] <- mvtnorm::rmvt(
        n = n_per_comp[k],
        sigma = mixture_params$Sigma_list[[k]],
        df = df,
        delta = mixture_params$mus_list[[k]],
        type = "shifted"
      )
    }
  }

  # Mixture log-density in unconstrained space
  log_density_unc_val <- log_density_mixture_unc(samps_unc, mixture_params)

  # Split into theta (unconstrained) and eta (free)
  theta_cols <- which(param_names %in% theta_names)
  eta_cols   <- which(!param_names %in% theta_names)

  theta_unc <- samps_unc[, theta_cols, drop = FALSE]
  eta_free  <- samps_unc[, eta_cols,   drop = FALSE]

  # Transform theta -> constrained
  trans_res <- transform_unconstrained_to_theta(theta_unc, theta_names,
                                                phi_min, phi_max)
  # Density in constrained space: log q(theta, eta) = log q_unc(u) - log|J|
  log_density <- log_density_unc_val - trans_res$log_jac

  # Expand seasonal states (sum-to-zero constraint)
  s_cols_free <- grep("^s\\d+$", colnames(eta_free), value = TRUE)
  if (length(s_cols_free) > 0) {
    last_s <- -rowSums(eta_free[, s_cols_free, drop = FALSE])
    eta_con <- cbind(eta_free, last_s)
    m_idx <- max(as.integer(sub("s", "", s_cols_free))) + 1
    colnames(eta_con)[ncol(eta_con)] <- paste0("s", m_idx)
  } else {
    eta_con <- eta_free
  }

  list(
    theta = trans_res$theta,
    eta = eta_con,
    theta_unc = theta_unc,
    eta_free = eta_free,
    log_density = log_density,
    log_density_unc = log_density_unc_val,
    log_jac = trans_res$log_jac
  )
}


##############################################################################
### AMIS Weight Calculation ###

# Compute AMIS weights using the balance heuristic over all historical proposals
# log_numerator: log(target_unc(u_i)) = log_target + log_jac for all particles
# all_samps_unc: matrix (N_total x d) of all unconstrained samples
# proposal_history: list of mixture_params, one per iteration
# n_per_iter: vector of sample sizes per iteration
# Returns: list(w = normalized weights, log_w = unnormalized log-weights)
compute_amis_weights <- function(log_numerator, all_samps_unc,
                                 proposal_history, n_per_iter) {
  N_total <- nrow(all_samps_unc)
  t_max <- length(proposal_history)

  # Balance-heuristic mixture weights: N_t / N_total
  log_mix_w <- log(n_per_iter) - log(N_total)

  # Evaluate each historical proposal at all particles
  log_q_matrix <- matrix(NA_real_, nrow = N_total, ncol = t_max)
  for (j in seq_len(t_max)) {
    log_q_matrix[, j] <- log_density_mixture_unc(
      all_samps_unc, proposal_history[[j]]
    ) + log_mix_w[j]
  }

  # Log-sum-exp across proposals: log(sum_t (N_t/N_total) * q_t(x_i))
  max_log_q <- apply(log_q_matrix, 1, max)
  log_denom <- max_log_q + log(rowSums(exp(log_q_matrix - max_log_q)))

  # Unnormalized log-weights
  log_w <- log_numerator - log_denom
  log_w[!is.finite(log_w)] <- -Inf

  # Self-normalize
  max_log_w <- max(log_w)
  w <- exp(log_w - max_log_w)
  w <- w / sum(w)

  list(w = w, log_w = log_w)
}


##############################################################################
### Weighted EM for MVT Mixture ###

# E-step computes responsibilities gamma_{i,k} = w_is[i] * r_{i,k}
# M-step updates pi_k, mu_k, Sigma_k with ridge regularization
# samps: matrix (N x d) of unconstrained samples
# w_is: self-normalized importance weight vector (sums to 1)
# mixture_params: current mixture parameters
# ridge_eps: ridge regularization added to diagonal of Sigma_k
# n_iter: number of EM iterations
em_update_mixture <- function(samps, w_is, mixture_params,
                              ridge_eps = 1e-4, n_iter = 5) {
  N <- nrow(samps)
  d <- ncol(samps)
  K <- mixture_params$K
  df <- mixture_params$df
  param_names <- colnames(samps)

  weights <- mixture_params$weights
  mus_list <- mixture_params$mus_list
  Sigma_list <- mixture_params$Sigma_list

  ridge_mat <- ridge_eps * diag(d)

  for (em_it in seq_len(n_iter)) {
    # ---- E-step: responsibilities ----
    log_resp <- matrix(NA_real_, nrow = N, ncol = K)
    for (k in seq_len(K)) {
      log_resp[, k] <- log(pmax(weights[k], 1e-300)) +
        mvtnorm::dmvt(samps,
                      delta = mus_list[[k]],
                      sigma = Sigma_list[[k]],
                      df = df, log = TRUE)
    }

    # Normalize across components (log-sum-exp per row)
    max_lr <- apply(log_resp, 1, max)
    resp <- exp(log_resp - max_lr)
    resp <- resp / rowSums(resp)

    # Incorporate importance weights: gamma_{i,k} = w_i * r_{i,k}
    gamma <- resp * w_is        # N x K
    gamma_sum <- colSums(gamma) # K-vector
    total_gamma <- sum(gamma_sum)

    # ---- M-step: mixing weights ----
    weights <- pmax(gamma_sum / total_gamma, 1e-6)
    weights <- weights / sum(weights)

    # ---- M-step: means and covariances ----
    for (k in seq_len(K)) {
      if (gamma_sum[k] < 1e-10) next

      gk <- gamma[, k]
      gk_sum <- gamma_sum[k]

      # Weighted mean
      mus_list[[k]] <- colSums(gk * samps) / gk_sum

      # Weighted covariance + ridge
      centered <- sweep(samps, 2, mus_list[[k]], "-")
      Sigma_list[[k]] <- crossprod(centered * sqrt(gk)) / gk_sum + ridge_mat

      # Enforce symmetry
      Sigma_list[[k]] <- 0.5 * (Sigma_list[[k]] + t(Sigma_list[[k]]))
    }
  }

  for (k in seq_len(K)) {
    names(mus_list[[k]]) <- param_names
    rownames(Sigma_list[[k]]) <- colnames(Sigma_list[[k]]) <- param_names
  }

  list(
    K = K, weights = weights,
    mus_list = mus_list, Sigma_list = Sigma_list,
    df = df
  )
}


##############################################################################
### AMIS Main Loop ###

adaptive_mis <- function(y, model_components, ctrl,
                         return_pointwise = FALSE,
                         use_nmig = FALSE) {

  # ---- Extract control parameters ----
  N_iter_max <- ctrl$N_iter_max
  N_draw     <- ctrl$N_draw
  N_draw_max <- ctrl$N_draw_max
  N_final    <- ctrl$N_final
  nu0        <- ctrl$nu0
  psi0       <- ctrl$psi0
  phi_min    <- ctrl$phi_min
  phi_max    <- ctrl$phi_max
  min_ess    <- ctrl$min_ess
  eta_df     <- ctrl$eta_df
  eta_df_incr_per_iter <- ctrl$eta_df_incr_per_iter
  N_draw_mult         <- ctrl$N_draw_mult
  first_iter_mult_N   <- ctrl$first_iter_mult_N
  verbose    <- ctrl$verbose
  v_spike    <- ctrl$v_spike
  v_slab     <- ctrl$v_slab
  w_nmig     <- ctrl$w_nmig
  K_mix        <- ctrl$K_mix
  ridge_eps    <- ctrl$ridge_eps
  em_iter      <- ctrl$em_iter
  jitter_scale <- ctrl$jitter_scale
  c_inflate_eta <- ctrl$c_inflate_eta

  L <- length(y)
  m <- stats::frequency(y)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha",
                    if (trend) c("beta", if (damped) "phi"),
                    if (seas) "gamma")

  # ---- Initialize ----
  base_params <- init_joint_params(y, model_components, theta_names, eta_df)
  mixture <- init_mixture_proposal(base_params, K = K_mix,
                                   jitter_scale = jitter_scale)

  eta_names <- setdiff(names(base_params$mus), theta_names)
  prior_eta_params <- list(
    mus   = base_params$mus[eta_names],
    Sigma = base_params$Sigma[eta_names, eta_names, drop = FALSE] * c_inflate_eta,
    df    = base_params$df
  )

  dummy_theta <- matrix(0, nrow = 1, ncol = length(theta_names))
  colnames(dummy_theta) <- theta_names
  log_prior_theta_const <- log_prior_theta_uniform(dummy_theta,
                                                   phi_min, phi_max)[1]

  # ---- AMIS accumulators ----
  all_samps_unc_list <- list()
  all_log_num_list   <- list()
  all_theta_list     <- list()
  all_eta_list       <- list()
  all_states_list    <- list()
  all_rss_list       <- list()
  proposal_history   <- list()
  n_per_iter         <- integer(0)

  timing <- list(draw = 0, refit = 0, weight = 0, update = 0, post = 0)
  amis_res <- NULL

  for (iter in seq_len(N_iter_max)) {

    # ---- Draw from current mixture ----
    t0 <- proc.time()[3]
    draws <- draw_from_mixture_proposal(N_draw, mixture, theta_names,
                                        phi_min, phi_max)
    timing$draw <- timing$draw + (proc.time()[3] - t0)

    # ---- Evaluate model (RSS + states) ----
    t0 <- proc.time()[3]
    refit <- RSS_vect_arma(
      yR = as.numeric(y),
      trend = trend, seas = seas, damped = damped, m = m,
      init_statesR = draws$eta,
      paramsR = draws$theta,
      return_residuals = FALSE
    )
    rss <- c(refit$RSS)
    timing$refit <- timing$refit + (proc.time()[3] - t0)

    # ---- Compute log target ----
    t0 <- proc.time()[3]
    log_lik <- -(nu0 + L) / 2 * log(psi0 + rss)
    log_prior_eta <- mvtnorm::dmvt(
      x     = draws$eta_free,
      delta = prior_eta_params$mus,
      sigma = prior_eta_params$Sigma,
      df    = prior_eta_params$df,
      log   = TRUE
    )
    log_target <- log_lik + log_prior_eta + log_prior_theta_const

    if (use_nmig) {
      nmig_cols <- grep("^(b|s\\d+)$", colnames(draws$eta_free), value = TRUE)
      if (length(nmig_cols) > 0) {
        nmig_vals <- draws$eta_free[, nmig_cols, drop = FALSE]
        log_nmig <- rowSums(vapply(
          seq_len(ncol(nmig_vals)),
          function(j) log_prior_nmig(nmig_vals[, j], v_spike, v_slab, w_nmig),
          numeric(nrow(nmig_vals))
        ))
        log_target <- log_target + log_nmig
      }
    }

    # ---- Accumulate iteration data ----
    all_samps_unc_list[[iter]] <- cbind(draws$theta_unc, draws$eta_free)
    all_log_num_list[[iter]]   <- log_target + draws$log_jac
    all_theta_list[[iter]]     <- draws$theta
    all_eta_list[[iter]]       <- draws$eta
    all_states_list[[iter]]    <- refit$states
    all_rss_list[[iter]]       <- rss
    proposal_history[[iter]]   <- mixture
    n_per_iter <- c(n_per_iter, N_draw)

    # ---- AMIS weight computation (all historical particles) ----
    all_samps_unc_mat <- do.call(rbind, all_samps_unc_list)
    all_log_num_vec   <- unlist(all_log_num_list)

    amis_res <- compute_amis_weights(all_log_num_vec, all_samps_unc_mat,
                                     proposal_history, n_per_iter)
    w   <- amis_res$w
    ess <- 1 / sum(w^2)
    timing$weight <- timing$weight + (proc.time()[3] - t0)

    if (verbose >= 2) {
      cat(sprintf("\n\nAdaptive Multiple IS - iter %d\n", iter))
      cat(sprintf("\n ESS = %.1f  (N_total = %d)\n", ess, sum(n_per_iter)))
    }

    if (ess >= min_ess) break
    if (iter == N_iter_max) {
      warning("maximum number of AMIS iterations reached")
      break
    }

    # ---- Update mixture via weighted EM ----
    t0 <- proc.time()[3]
    mixture <- em_update_mixture(all_samps_unc_mat, w, mixture,
                                 ridge_eps = ridge_eps, n_iter = em_iter)
    mixture$df <- mixture$df + eta_df_incr_per_iter

    if (iter >= first_iter_mult_N) {
      N_draw <- min(as.integer(N_draw * N_draw_mult), N_draw_max)
    }
    timing$update <- timing$update + (proc.time()[3] - t0)
  }

  # ---- Post-processing ----
  t0 <- proc.time()[3]
  N_total <- sum(n_per_iter)
  all_theta_mat  <- do.call(rbind, all_theta_list)
  all_eta_mat    <- do.call(rbind, all_eta_list)
  all_states_mat <- do.call(rbind, all_states_list)
  all_rss_vec    <- unlist(all_rss_list)
  colnames(all_states_mat) <- colnames(all_eta_mat)

  res_idx <- sample(N_total, size = N_final, replace = TRUE, prob = w)
  thetas  <- all_theta_mat[res_idx, , drop = FALSE]
  etas    <- all_eta_mat[res_idx, , drop = FALSE]
  states  <- all_states_mat[res_idx, , drop = FALSE]
  sigma2s <- (psi0 + all_rss_vec[res_idx]) / stats::rchisq(N_final, df = nu0 + L)

  log_evidence <- max(amis_res$log_w) +
    log(mean(exp(amis_res$log_w - max(amis_res$log_w))))

  log_lik_pointwise <- NULL
  if (return_pointwise) {
    refit_final <- RSS_vect_arma(
      yR = as.numeric(y),
      trend = trend, seas = seas, damped = damped, m = m,
      init_statesR = etas,
      paramsR = thetas,
      return_residuals = TRUE
    )
    E <- refit_final$residuals
    sd_mat <- matrix(sqrt(sigma2s), nrow = nrow(E), ncol = ncol(E),
                     byrow = FALSE)
    log_lik_pointwise <- stats::dnorm(E, mean = 0, sd = sd_mat, log = TRUE)
  }
  timing$post <- timing$post + (proc.time()[3] - t0)

  list(
    thetas = thetas,
    etas = etas,
    states = states,
    sigma2s = sigma2s,
    ess = ess,
    n_iter = iter,
    prop_params = mixture,
    log_evidence = log_evidence,
    log_lik_pointwise = log_lik_pointwise,
    timing = timing
  )
}
