################################################################################
# ADAPTIVE IMPORTANCE SAMPLING (AIS)
#
# Samples the smoothing parameters only: the initial states and \sigma^2 are
# integrated analytically (conjugate prior). The proposal is a Student-t at the
# posterior mode, adapted while the ESS is below min_ess.

# Fit of one model by AIS
adaptive_is_rb <- function(y, model_components, ctrl,
                           return_pointwise = FALSE) {
  N_iter_max  <- ctrl$N_iter_max
  N_draw_raw  <- ctrl$N_draw
  N_final_raw <- ctrl$N_final
  min_ess_raw <- ctrl$min_ess   # NULL: N_draw / 4
  is_df       <- ctrl$is_df
  is_scale    <- ctrl$is_scale
  lr          <- ctrl$lr
  nu0         <- ctrl$nu0
  verbose     <- ctrl$verbose

  L <- length(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha", if (trend) c("beta", if (damped) "phi"), if (seas) "gamma")
  n_theta <- length(theta_names)

  # Values for d = n_theta
  N_draw  <- resolve_by_d(N_draw_raw,  n_theta)
  N_final <- resolve_by_d(N_final_raw, n_theta)
  min_ess <- if (is.null(min_ess_raw)) N_draw / 4 else resolve_by_d(min_ess_raw, n_theta)

  prior <- init_rb_prior(y, model_components, theta_names, ctrl)
  log_g_rb <- make_log_g_rb(y, model_components, theta_names, ctrl, prior)

  # AIS, with the mode search starting at the best point of a prior scan
  t0 <- proc.time()[3]
  scan <- prior_scan(log_g_rb, theta_names, ctrl$n_scan)
  t_scan <- proc.time()[3] - t0
  ais <- adaptive_importance_sampling(log_g_rb, scan$Z[1, ],
                                      n_draw = N_draw, min_ess = min_ess,
                                      df = is_df, scale = is_scale,
                                      n_iter_max = N_iter_max, lr = lr,
                                      verbose = verbose)
  timing <- c(list(scan = t_scan), ais$timing, list(post = 0))

  prop_params <- ais$proposal
  names(prop_params$mus) <- theta_names
  rownames(prop_params$Sigma) <- colnames(prop_params$Sigma) <- theta_names

  # Target ESS not reached: the model gets zero weight
  if (ais$ess < min_ess) {
    warning(sprintf(paste0(
      "adaptive_is_rb: ESS = %.0f below min_ess = %.0f after %d iterations; ",
      "the model gets zero weight"), ais$ess, min_ess, ais$n_iter))
    return(list(
      thetas = NULL,
      etas = NULL,
      states = NULL,
      sigma2s = NULL,
      ess = ais$ess,
      n_iter = ais$n_iter,
      prop_params = prop_params,
      log_evidence = -Inf,
      log_lik_pointwise = if (return_pointwise) matrix(-1e300, nrow = N_final, ncol = L) else NULL,  # not -Inf: NaN in logsumexp
      timing = timing
    ))
  }

  # Posterior draws, from the draws of all iterations
  t0 <- proc.time()[3]
  theta_particles <- do.call(rbind, lapply(ais$draw_evals, `[[`, "theta"))
  ml_res <- bind_ml_res(lapply(ais$draw_evals, `[[`, "ml_res"))
  post <- draw_rb_posterior(y, model_components, theta_particles, ais$w, ml_res,
                            N_final = N_final, nu0 = nu0,
                            return_pointwise = return_pointwise)
  timing$post <- proc.time()[3] - t0

  list(
    thetas = post$thetas,
    etas = post$etas,
    states = post$states,
    sigma2s = post$sigma2s,
    ess = ais$ess,
    n_iter = ais$n_iter,
    prop_params = prop_params,
    log_evidence = ais$log_evidence,
    log_lik_pointwise = post$log_lik_pointwise,
    timing = timing
  )
}


################################################################################
# PRIOR AND POSTERIOR DRAWS (shared with quadrature_rb())

# Prior of the initial states and \sigma^2, and log-density of the uniform
# prior of the smoothing parameters. eta0 and V0 are in the C++ order.
init_rb_prior <- function(y, model_components, theta_names, ctrl) {
  psi0          <- ctrl$psi0
  phi_min       <- ctrl$phi_min
  phi_max       <- ctrl$phi_max
  c_inflate_eta <- ctrl$c_inflate_eta

  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  n_theta <- length(theta_names)

  eta_init <- init_eta_params(y, model_components)
  eta_names_free <- names(eta_init$mus)  # l, [b,] [s1, ..., s_{m-1}]

  # eta also includes the last seasonal state: (l0, [b0,] [s1, ..., s_m])
  n_eta <- 1 + (if (trend) 1L else 0L) + (if (seas) m else 0L)

  # Prior mean of s_m: -(s1 + ... + s_{m-1})
  eta0_free <- eta_init$mus
  if (seas && m > 1) {
    s_free_names <- grep("^s\\d+$", names(eta0_free), value = TRUE)
    last_s <- -sum(eta0_free[s_free_names])
    eta0_r_order <- c(eta0_free, last_s)
    names(eta0_r_order)[length(eta0_r_order)] <- paste0("s", m)
  } else {
    eta0_r_order <- eta0_free
  }

  # eta | \sigma^2 ~ N(eta0, \sigma^2 V0), with V0 = Sigma / E[\sigma^2]: prior
  # covariance Sigma at the prior mean of \sigma^2
  Sigma_heuristic_free <- eta_init$Sigma * c_inflate_eta
  prior_mean_sigma2 <- psi0 / (ctrl$nu0 - 2)

  # Independent priors on all m seasonal states, without the sum-to-zero
  # constraint: the prior fixes the shift (l0 + c, s - c) that the data cannot
  # identify, and keeps V0 invertible
  if (seas && m > 1) {
    n_free <- length(eta0_free)
    n_full <- n_eta
    s_indices_free <- grep("^s\\d+$", names(eta0_free))
    Sigma_full <- matrix(0, n_full, n_full)
    Sigma_full[1:n_free, 1:n_free] <- Sigma_heuristic_free
    # s_m: same prior variance as the other seasonal states
    Sigma_full[n_full, n_full] <- mean(diag(Sigma_heuristic_free)[s_indices_free])
  } else {
    Sigma_full <- Sigma_heuristic_free
  }

  # C++ order of the seasonal states: (s_m, ..., s1)
  if (seas && m > 1) {
    s_offset_r <- (1 + (if (trend) 1L else 0L))  # number of non-seasonal states
    non_s_idx <- seq_len(s_offset_r)
    s_idx_r <- (s_offset_r + 1):(s_offset_r + m)  # s1, ..., s_m
    s_idx_cpp <- rev(s_idx_r)                     # s_m, ..., s1
    reorder <- c(non_s_idx, s_idx_cpp)

    eta0_cpp <- as.numeric(eta0_r_order[reorder])
    V0 <- Sigma_full[reorder, reorder, drop = FALSE] / prior_mean_sigma2
  } else {
    eta0_cpp <- as.numeric(eta0_r_order)
    V0 <- Sigma_full / prior_mean_sigma2
  }

  # Uniform prior of the smoothing parameters: constant log-density
  dummy_theta <- matrix(0, nrow = 1, ncol = n_theta)
  colnames(dummy_theta) <- theta_names
  log_prior_theta_const <- log_prior_theta_uniform(dummy_theta, phi_min, phi_max)[1]

  list(
    eta0 = eta0_cpp,
    V0 = V0,
    psi0 = psi0,
    n_eta = n_eta,
    log_prior_theta_const = log_prior_theta_const
  )
}

# Posterior draws: resampling of the weighted particles, then \sigma^2 and eta
# from their conditional posteriors. ml_res: output of marginal_likelihood_rb()
# at the particles.
draw_rb_posterior <- function(y, model_components, theta_particles, w, ml_res,
                              N_final, nu0, return_pointwise = FALSE) {
  L <- length(y)
  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  n_eta <- 1 + (if (trend) 1L else 0L) + (if (seas) m else 0L)
  N_particles <- nrow(theta_particles)

  res_idx <- sample(N_particles, size = N_final, replace = TRUE, prob = w)
  thetas <- theta_particles[res_idx, , drop = FALSE]

  posterior_scale <- as.numeric(ml_res$posterior_scale)
  nu_n <- nu0 + L

  sigma2s <- posterior_scale[res_idx] / stats::rchisq(N_final, df = nu_n)

  # eta ~ N(mu_n, \sigma^2 V_n), with V_n = t(Rn) Rn; all draws at once
  # (column j of Z: standard normals of draw j)
  Z <- matrix(stats::rnorm(n_eta * N_final), nrow = n_eta, ncol = N_final)
  Rn_draws <- ml_res$Rn[, , res_idx, drop = FALSE]
  Z_draws  <- array(Z[, rep(seq_len(N_final), each = n_eta)], dim = c(n_eta, n_eta, N_final))
  etas <- t(ml_res$mu_n[, res_idx, drop = FALSE] +
              sweep(colSums(Rn_draws * Z_draws), 2, sqrt(sigma2s), "*"))

  # C++ order: l, [b,] s_m, ..., s1
  eta_col_names_cpp <- c("l", if (trend) "b",
                         if (seas) paste0("s", rev(seq_len(m))))
  colnames(etas) <- eta_col_names_cpp

  # R order: l, [b,] s1, ..., s_m
  eta_col_names_r <- c("l", if (trend) "b",
                       if (seas) paste0("s", seq_len(m)))
  etas <- etas[, eta_col_names_r, drop = FALSE]

  # Final states
  refit_final <- RSS_vect_arma(
    yR = as.numeric(y),
    trend = trend,
    seas = seas,
    damped = damped,
    m = m,
    init_statesR = etas,
    paramsR = thetas,
    return_residuals = return_pointwise
  )
  states <- refit_final$states
  colnames(states) <- eta_col_names_r

  log_lik_pointwise <- NULL
  if (return_pointwise) {
    E <- refit_final$residuals
    sd_mat <- matrix(sqrt(sigma2s), nrow = nrow(E), ncol = ncol(E), byrow = FALSE)
    log_lik_pointwise <- stats::dnorm(E, mean = 0, sd = sd_mat, log = TRUE)
  }

  list(
    thetas = thetas,
    etas = etas,
    states = states,
    sigma2s = sigma2s,
    log_lik_pointwise = log_lik_pointwise
  )
}


################################################################################
# AIS HELPERS

# AIS of g = exp(log g) on R^d (log_g_fn, z_start as in adaptive_gh_quadrature()).
# Student-t proposal at the mode of log g, with scale matrix scale * H^{-1};
# updated while ESS < min_ess, adding n_draw draws per iteration. The draws of
# all iterations are pooled and weighted with the mixture of the proposals.
adaptive_importance_sampling <- function(log_g_fn, z_start, n_draw, min_ess,
                                         df = 5, scale = 4,
                                         n_iter_max = 30, lr = 0.9,
                                         verbose = 0, ...) {
  lap <- laplace_mode(log_g_fn, z_start, ...)
  proposal <- list(mus = lap$zhat, Sigma = scale * tcrossprod(lap$Lmat), df = df)
  timing <- c(lap$timing, list(draws = 0, update = 0))

  proposals <- list()
  draw_evals <- list()
  Z <- NULL
  log_g <- NULL
  for (iter in seq_len(n_iter_max)) {
    # New draws, evaluated in one batch
    t0 <- proc.time()[3]
    Z_new <- draw_t_rqmc(n_draw, proposal$mus, proposal$Sigma, proposal$df)
    draw_eval <- log_g_fn(Z_new)
    timing$draws <- timing$draws + (proc.time()[3] - t0)
    proposals[[iter]]  <- c(proposal, list(n = nrow(Z_new)))
    draw_evals[[iter]] <- draw_eval
    Z <- rbind(Z, Z_new)
    log_g <- c(log_g, draw_eval$log_g)

    # Weights of the pooled draws
    log_w <- log_g - log_mix_density(Z, proposals)
    log_w[!is.finite(log_w)] <- -Inf
    lw_max <- max(log_w)
    if (!is.finite(lw_max)) {
      w <- NULL
      log_evidence <- -Inf
      ess <- 0
      break
    }
    w <- exp(log_w - lw_max)
    log_evidence <- lw_max + log(mean(w))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)

    if (isTRUE(verbose >= 2)) {
      cat(sprintf("\n\nRao-Blackwellized AIS - iter %d\n", iter))
      cat(sprintf("\n ESS = %.1f\n", ess))
    }
    if (ess >= min_ess || iter == n_iter_max) break

    # Proposal update
    t0 <- proc.time()[3]
    proposal <- update_theta_only_proposal(Z, w, proposal, lr = lr)
    timing$update <- timing$update + (proc.time()[3] - t0)
  }

  list(
    log_evidence = log_evidence,
    w = w,
    Z = Z,
    draw_evals = draw_evals,
    ess = ess,
    n_iter = iter,
    proposal = proposal,
    zhat = lap$zhat,
    timing = timing
  )
}

# n draws from a multivariate Student-t, from randomised Sobol points
draw_t_rqmc <- function(n, mu, Sigma, df) {
  d <- length(mu)
  u <- matrix(qrng::sobol(n, d = d + 1, randomize = "digital.shift"), ncol = d + 1)
  eps  <- stats::qnorm(u[, seq_len(d), drop = FALSE])
  chi2 <- stats::qchisq(u[, d + 1], df = df)
  devs <- (eps %*% chol(Sigma)) / sqrt(chi2 / df)
  # Drop draws from points exactly on 0 or 1
  devs <- devs[rowSums(!is.finite(devs)) == 0, , drop = FALSE]
  sweep(devs, 2, mu, "+")
}

# Log-density of the mixture of the proposals at the rows of Z (weights: share
# of the draws)
log_mix_density <- function(Z, proposals) {
  n_tot <- sum(vapply(proposals, `[[`, numeric(1), "n"))
  log_q <- vapply(proposals, function(p) {
    log(p$n / n_tot) + ldmvt_chol(sweep(Z, 2, p$mus, "-"), chol(p$Sigma), p$df)
  }, numeric(nrow(Z)))
  log_q <- matrix(log_q, nrow = nrow(Z))
  q_max <- apply(log_q, 1, max)
  q_max + log(rowSums(exp(log_q - q_max)))
}

# Binds the outputs of marginal_likelihood_rb() for several batches of particles
bind_ml_res <- function(ml_list) {
  if (length(ml_list) == 1L) return(ml_list[[1]])
  n_eta <- nrow(ml_list[[1]]$mu_n)
  n_tot <- sum(vapply(ml_list, function(r) length(r$posterior_scale), numeric(1)))
  list(
    log_marginal_lik = unlist(lapply(ml_list, function(r) as.numeric(r$log_marginal_lik))),
    Rn               = array(unlist(lapply(ml_list, `[[`, "Rn")), dim = c(n_eta, n_eta, n_tot)),
    mu_n             = do.call(cbind, lapply(ml_list, `[[`, "mu_n")),
    posterior_scale  = unlist(lapply(ml_list, function(r) as.numeric(r$posterior_scale)))
  )
}

# Proposal update by weighted moment matching
update_theta_only_proposal <- function(theta_unc, w, prev_params,
                                       lr = 0.9, min_var = 1e-6,
                                       lambda_shr = 0.1) {
  w <- w / sum(w)
  ess <- 1 / sum(w^2)

  new_mus <- colSums(w * theta_unc)

  centered <- sweep(theta_unc, 2, new_mus, "-")
  Sigma_new <- crossprod(centered * sqrt(w))

  diag_Sigma <- pmax(diag(Sigma_new), min_var)
  Sigma_new <- (1 - lambda_shr) * Sigma_new
  diag(Sigma_new) <- diag_Sigma

  lamb <- min(lr, ess / (100 + ess))

  out_mus <- lamb * new_mus + (1 - lamb) * prev_params$mus
  out_Sigma <- lamb * Sigma_new + (1 - lamb) * prev_params$Sigma

  list(mus = out_mus, Sigma = out_Sigma, df = prev_params$df)
}
