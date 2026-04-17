##############################################################################
### Adaptive Importance Sampling ###

adaptive_is <- function(y, model_components, ctrl,
                        return_pointwise = FALSE,
                        use_nmig = FALSE) {
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
  lr         <- ctrl$lr
  c_inflate_eta       <- ctrl$c_inflate_eta
  N_draw_mult         <- ctrl$N_draw_mult
  first_iter_mult_N   <- ctrl$first_iter_mult_N
  factor_inflate_Sigma <- ctrl$factor_inflate_Sigma
  verbose    <- ctrl$verbose
  v_spike    <- ctrl$v_spike
  v_slab     <- ctrl$v_slab
  w_nmig     <- ctrl$w_nmig

  L <- length(y)
  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha", if (trend) c("beta", if (damped) "phi"), if (seas) "gamma")

  prop_params <- init_joint_params(y, model_components, theta_names, eta_df)

  eta_names <- setdiff(colnames(prop_params$Sigma), theta_names)
  prior_eta_params <- list(
    mus = prop_params$mus[eta_names],
    Sigma = prop_params$Sigma[eta_names, eta_names, drop = FALSE] * c_inflate_eta,
    df = prop_params$df
  )

  dummy_theta <- matrix(0, nrow = 1, ncol = length(theta_names))
  colnames(dummy_theta) <- theta_names
  log_prior_theta_const <- log_prior_theta_uniform(dummy_theta, phi_min, phi_max)[1]

  prev_ess <- 0
  timing <- list(draw = 0, refit = 0, weight = 0, update = 0, post = 0)
  do_time <- verbose >= 1

  # Pre-compute once: avoid repeated as.numeric(y) and Cholesky inside the loop
  y_vec        <- as.numeric(y)
  R_chol       <- chol(prop_params$Sigma)               # proposal Cholesky
  R_prior_eta  <- chol(prior_eta_params$Sigma)          # prior-eta Cholesky (constant)
  d_eta        <- length(eta_names)
  df_eta       <- prior_eta_params$df
  mu_eta       <- prior_eta_params$mus
  log_det_R_pe <- sum(log(diag(R_prior_eta)))
  lc_eta       <- lgamma((df_eta + d_eta) / 2) - lgamma(df_eta / 2) -
                   (d_eta / 2) * log(df_eta * pi) - log_det_R_pe

  for (iter in seq_len(N_iter_max)) {
    if (do_time) t0 <- proc.time()[3]
    draws <- draw_from_joint_proposal(N_draw, prop_params, theta_names, phi_min, phi_max,
                                      chol_Sigma = R_chol)
    theta_samp <- draws$theta
    eta_samp <- draws$eta
    if (do_time) timing$draw <- timing$draw + (proc.time()[3] - t0)

    if (do_time) t0 <- proc.time()[3]
    refit <- RSS_vect_arma(
      yR = y_vec,
      trend = trend,
      seas = seas,
      damped = damped,
      m = m,
      init_statesR = eta_samp,
      paramsR = theta_samp,
      return_residuals = return_pointwise
    )
    rss <- c(refit$RSS)
    if (do_time) timing$refit <- timing$refit + (proc.time()[3] - t0)

    if (do_time) t0 <- proc.time()[3]
    log_lik <- -(nu0 + L) / 2 * log(psi0 + rss)
    # Prior log-density for eta: MVT(mu_eta, Sigma_eta, df_eta)
    # lc_eta and R_prior_eta are constant across iterations — computed once above.
    devs_eta      <- sweep(draws$eta_free, 2, mu_eta, "-")
    z_eta         <- forwardsolve(t(R_prior_eta), t(devs_eta))
    mahal_eta     <- colSums(z_eta^2)
    log_prior_eta <- lc_eta - ((df_eta + d_eta) / 2) * log(1 + mahal_eta / df_eta)
    log_target <- log_lik + log_prior_eta + log_prior_theta_const

    # --- NMIG spike-and-slab penalty on trend / seasonal initial states ---
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

    log_w <- log_target - draws$log_density
    w <- exp(log_w - max(log_w))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)
    if (do_time) timing$weight <- timing$weight + (proc.time()[3] - t0)

    if (verbose >= 2) cat(sprintf("\n\nAdaptive Importance Sampling - iter %d\n", iter))
    if (verbose >= 2) cat(sprintf("\n ESS = %.1f\n", ess))

    if (ess >= min_ess) break
    if (iter == N_iter_max) {
      warning("maximum number of AIS iterations reached")
      break
    }

    if (do_time) t0 <- proc.time()[3]
    prop_params <- update_joint_proposal(
      theta_unc = draws$theta_unc,
      eta_free = draws$eta_free,
      w = w,
      prev_params = prop_params,
      lr = lr
    )

    prop_params$df <- prop_params$df + eta_df_incr_per_iter
    if (iter >= first_iter_mult_N) {
      N_draw <- min(as.integer(N_draw * N_draw_mult), N_draw_max)
    }
    if (iter > 1 && ess < 0.8 * prev_ess) {
      prop_params$Sigma <- prop_params$Sigma * factor_inflate_Sigma
    }
    R_chol <- chol(prop_params$Sigma)   # recompute once after Sigma update
    if (do_time) timing$update <- timing$update + (proc.time()[3] - t0)
    prev_ess <- ess
  }

  if (do_time) t0 <- proc.time()[3]
  res_idx <- sample(N_draw, size = N_final, replace = TRUE, prob = w)
  thetas <- theta_samp[res_idx, , drop = FALSE]
  etas <- eta_samp[res_idx, , drop = FALSE]
  states <- refit$states[res_idx, , drop = FALSE]
  colnames(states) <- colnames(etas)
  sigma2s <- (psi0 + rss[res_idx]) / stats::rchisq(N_final, df = nu0 + L)

  log_evidence <- max(log_w) + log(mean(exp(log_w - max(log_w))))

  log_lik_pointwise <- NULL
  if (return_pointwise) {
    E <- refit$residuals[res_idx, , drop = FALSE]
    sd_mat <- matrix(sqrt(sigma2s), nrow = nrow(E), ncol = ncol(E), byrow = FALSE)
    log_lik_pointwise <- stats::dnorm(E, mean = 0, sd = sd_mat, log = TRUE)
  }
  if (do_time) timing$post <- timing$post + (proc.time()[3] - t0)

  list(
    thetas = thetas,
    etas = etas,
    states = states,
    sigma2s = sigma2s,
    ess = ess,
    n_iter = iter,
    prop_params = prop_params,
    log_evidence = log_evidence,
    log_lik_pointwise = log_lik_pointwise,
    timing = timing
  )
}

##############################################################################
### Stacking Optimization ###

compute_stacking_weights <- function(log_lik_list) {
  K <- length(log_lik_list)
  if (K == 0) stop("log_lik_list is empty")
  if (is.null(log_lik_list[[1]]) || !is.matrix(log_lik_list[[1]])) {
    stop("Each element of log_lik_list must be a matrix (samples x time)")
  }
  L <- ncol(log_lik_list[[1]])

  lpd_model <- matrix(nrow = L, ncol = K)
  for (k in seq_len(K)) {
    ll_mat <- log_lik_list[[k]]
    if (!is.matrix(ll_mat) || ncol(ll_mat) != L) {
      stop("All log_lik_list elements must be matrices with same number of columns")
    }
    S <- nrow(ll_mat)
    if (S <= 0) stop("Each log-likelihood matrix must have at least one row")
    max_ll <- apply(ll_mat, 2, max)
    sum_exp <- colSums(exp(ll_mat - matrix(max_ll, nrow = S, ncol = L, byrow = TRUE)))
    lpd_model[, k] <- max_ll + log(sum_exp) - log(S)
  }

  exp_lpd <- exp(lpd_model)

  obj_fun <- function(par) {
    w <- exp(par) / sum(exp(par))
    mix_dens <- as.vector(exp_lpd %*% w)
    -sum(log(pmax(mix_dens, 1e-300)))
  }

  opt <- stats::optim(rep(0, K), obj_fun, method = "BFGS")
  w_final <- exp(opt$par) / sum(exp(opt$par))

  w_final[w_final < 1e-3] <- 0
  if (sum(w_final) <= 0 || any(!is.finite(w_final))) {
    rep(1 / K, K)
  } else {
    w_final / sum(w_final)
  }
}

##############################################################################
### Model Fitting Wrapper ###

fit_bets_models <- function(y,
                            model_components,
                            ctrl,
                            method = c("bma", "stacking", "nmig"),
                            sampler = c("ais", "amis"),
                            rao_blackwellize_eta = FALSE) {
  method <- match.arg(method)
  sampler <- match.arg(sampler)

  verbose <- ctrl$verbose
  psi0    <- ctrl$psi0

  m <- stats::frequency(y)

  # --- NMIG: override model_components to a single super-model ------------
  use_nmig <- (method == "nmig")
  if (use_nmig) {
    if (m > 1) {
      model_components <- list(c("A", "A", "A", "TRUE"))   # ETS(A,Ad,A)
    } else {
      model_components <- list(c("A", "A", "N", "TRUE"))   # ETS(A,Ad,N)
    }
  }

  n_models <- length(model_components)
  need_pointwise <- (method == "stacking")

  if (is.null(psi0)) {
    mse_naive <- mean(diff(y, lag = 1)^2)
    psi0 <- if (m > 1) 0.5 * (mean(diff(y, lag = m)^2) + mse_naive) else mse_naive
    ctrl$psi0 <- psi0
  }

  results_list <- vector("list", n_models)
  log_marginal_liks <- rep(NA_real_, n_models)
  log_lik_list <- if (need_pointwise) vector("list", n_models) else NULL
  fit_time_per_model <- numeric(n_models)

  for (i in seq_along(model_components)) {
    if (verbose >= 2) cat(sprintf("\nFitting model %d of %d\n", i, n_models))
    t0 <- proc.time()[3]
    if (rao_blackwellize_eta) {
      sampler_fn <- adaptive_is_rb
    } else {
      sampler_fn <- if (sampler == "amis") adaptive_mis else adaptive_is
    }
    res_i <- sampler_fn(
      y,
      model_components[[i]],
      ctrl = ctrl,
      return_pointwise = need_pointwise,
      use_nmig = use_nmig
    )
    fit_time_per_model[i] <- proc.time()[3] - t0
    results_list[[i]] <- res_i
    results_list[[i]]$model_components <- model_components[[i]]
    log_marginal_liks[i] <- res_i$log_evidence
    if (need_pointwise) log_lik_list[[i]] <- res_i$log_lik_pointwise
  }

  t0 <- proc.time()[3]
  if (use_nmig) {
    model_weights <- 1.0
    if (verbose >= 1) {
      cat("\nNMIG Spike-and-Slab (single super-model):\n")
      cat(sprintf("  %-5s: %.3f\n", ets_label(model_components[[1]]), 1.0))
    }
  } else if (need_pointwise) {
    model_weights <- compute_stacking_weights(log_lik_list)
    if (verbose >= 1) cat("\nStacking Weights:\n")
  } else {
    prior_models <- ctrl$prior_models
    if (is.null(prior_models)) prior_models <- rep(1 / n_models, n_models)
    log_post_unnorm <- log(prior_models) + log_marginal_liks
    model_weights <- exp(log_post_unnorm - max(log_post_unnorm))
    model_weights <- model_weights / sum(model_weights)
    if (verbose >= 1) cat("\nBMA Weights:\n")
  }
  weight_time <- proc.time()[3] - t0

  if (verbose >= 1 && !use_nmig) {
    labels <- vapply(model_components, ets_label, character(1))
    for (i in seq_len(n_models)) cat(sprintf("  %-5s: %.3f\n", labels[i], model_weights[i]))
  }

  list(
    results = results_list,
    log_marginal_liks = log_marginal_liks,
    model_weights = model_weights,
    timing = list(
      fit_time_per_model = fit_time_per_model,
      fit_time_total = sum(fit_time_per_model),
      weight_time = weight_time,
      total_time = sum(fit_time_per_model) + weight_time
    )
  )
}

##############################################################################
### Trajectory Simulation ###

simulate_future_trajectories <- function(bets_fit, h = 10, n_traj = 1000) {
  n_models <- length(bets_fit$results)
  model_weights <- bets_fit$model_weights

  traj_list <- vector("list", n_models)
  n_traj_list <- round(n_traj * model_weights)
  n_traj_list[which.max(n_traj_list)] <- n_traj - sum(n_traj_list) + max(n_traj_list)

  for (i in seq_along(bets_fit$results)) {
    if (n_traj_list[i] > 0) {
      idxs <- sample(nrow(bets_fit$results[[i]]$thetas), size = n_traj_list[i], replace = TRUE)
      res_i <- bets_fit$results[[i]]
      traj_list[[i]] <- ets_future_traj(
        model_components = res_i$model_components,
        states = res_i$states[idxs, , drop = FALSE],
        params = res_i$thetas[idxs, , drop = FALSE],
        sigma2s = res_i$sigma2s[idxs],
        h = h
      )
    }
  }

  do.call(rbind, traj_list)
}

##############################################################################
### Helpers ###

ets_label <- function(model_components) {
  d <- if (model_components[[4]] == "TRUE" && model_components[[2]] == "A") "d" else ""
  paste0(model_components[[1]], model_components[[2]], d, model_components[[3]])
}
