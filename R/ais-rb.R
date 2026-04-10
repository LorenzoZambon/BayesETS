##############################################################################
### Rao-Blackwellized Adaptive Importance Sampling ###
###
### Samples only theta, analytically integrating out eta (initial states)
### and sigma^2 using conjugate Normal-Inverse-Chi-Squared priors.
##############################################################################

adaptive_is_rb <- function(y, model_components, ctrl,
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

  L <- length(y)
  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha", if (trend) c("beta", if (damped) "phi"), if (seas) "gamma")
  n_theta <- length(theta_names)

  # ---- Set up eta prior ----
  # Use the existing heuristic for the initial states
  eta_init <- init_eta_params(y, model_components, eta_df = eta_df)
  eta_names_free <- names(eta_init$mus)  # l, [b,] [s1..s_{m-1}]

  # For the RB formulation, eta includes the last seasonal state too
  # eta = (l0, [b0,] [s1, ..., s_m])  -- 1-indexed; s_m is the sum-to-zero state
  n_eta <- 1 + (if (trend) 1L else 0L) + (if (seas) m else 0L)

  # Expand eta0 to include s_m = -sum(s1..s_{m-1})
  eta0_free <- eta_init$mus
  if (seas && m > 1) {
    s_free_names <- grep("^s\\d+$", names(eta0_free), value = TRUE)
    last_s <- -sum(eta0_free[s_free_names])
    eta0_r_order <- c(eta0_free, last_s)
    names(eta0_r_order)[length(eta0_r_order)] <- paste0("s", m)
  } else {
    eta0_r_order <- eta0_free
  }

  # Build full prior covariance matching the heuristic
  # E[sigma^2] = psi0/(nu0-2) = psi0 (since nu0=3)
  # So V0 = Sigma_heuristic / psi0
  Sigma_heuristic_free <- eta_init$Sigma * c_inflate_eta

  # Expand to include the constrained seasonal state.
  # The sum-to-zero constraint s_1 + ... + s_m = 0 is structural (enforced by
  # the design matrix X), so the prior does not need to encode it.  Encoding
  # it via off-diagonal coupling terms makes Sigma_full exactly rank-deficient
  # (the direction [0,...,0,1,...,1] over all seasonal states has zero
  # quadratic form), causing inv_sympd(V0) to fail.  Use an independent prior
  # for s_m with the same marginal scale as the other seasonal states instead.
  if (seas && m > 1) {
    n_free <- length(eta0_free)
    n_full <- n_eta
    s_indices_free <- grep("^s\\d+$", names(eta0_free))
    Sigma_full <- matrix(0, n_full, n_full)
    Sigma_full[1:n_free, 1:n_free] <- Sigma_heuristic_free
    # Independent prior for the constrained seasonal state s_m
    Sigma_full[n_full, n_full] <- mean(diag(Sigma_heuristic_free)[s_indices_free])
  } else {
    Sigma_full <- Sigma_heuristic_free
  }

  # build_design_and_c_batch maps eta[s_offset + k] directly to times t ≡ k (mod m).
  # Slot k=0 is used at t=0, m, 2m, … (oldest periodic factor = s_m in R naming).
  # Slot k=m-1 is used at t=m-1, 2m-1, … (most-recent factor = s_1 in R naming).
  # So the C++ ordering is: eta_cpp = (l, [b,] s_m, s_{m-1}, ..., s_1)
  # We must reverse the seasonal block of eta0 and V0 before passing to C++.
  if (seas && m > 1) {
    s_offset_r <- (1 + (if (trend) 1L else 0L))  # 1-based offset to first seasonal in R order
    non_s_idx <- seq_len(s_offset_r)
    s_idx_r <- (s_offset_r + 1):(s_offset_r + m)  # s1..s_m in R order
    s_idx_cpp <- rev(s_idx_r)  # s_m, s_{m-1}, ..., s1
    reorder <- c(non_s_idx, s_idx_cpp)

    eta0_cpp <- as.numeric(eta0_r_order[reorder])
    V0 <- Sigma_full[reorder, reorder, drop = FALSE] / psi0
  } else {
    eta0_cpp <- as.numeric(eta0_r_order)
    V0 <- Sigma_full / psi0
  }
  eta0 <- eta0_cpp

  # ---- Initialize theta-only proposal ----
  # Start from the joint init but extract theta portion only
  joint_init <- init_joint_params(y, model_components, theta_names, eta_df)
  theta_prop_params <- list(
    mus = joint_init$mus[theta_names],
    Sigma = joint_init$Sigma[theta_names, theta_names, drop = FALSE],
    df = joint_init$df
  )

  dummy_theta <- matrix(0, nrow = 1, ncol = n_theta)
  colnames(dummy_theta) <- theta_names
  log_prior_theta_const <- log_prior_theta_uniform(dummy_theta, phi_min, phi_max)[1]

  prev_ess <- 0
  timing <- list(draw = 0, design = 0, marglik = 0, weight = 0, update = 0, post = 0)

  for (iter in seq_len(N_iter_max)) {
    # ---- Step 2: Draw theta particles ----
    t0 <- proc.time()[3]
    draws <- draw_theta_only(N_draw, theta_prop_params, theta_names, phi_min, phi_max)
    timing$draw <- timing$draw + (proc.time()[3] - t0)

    # ---- Step 3: Build design matrices (C++) ----
    t0 <- proc.time()[3]
    design <- build_design_and_c_batch(
      yR = as.numeric(y),
      trend = trend,
      seas = seas,
      damped = damped,
      m = m,
      paramsR = draws$theta
    )
    timing$design <- timing$design + (proc.time()[3] - t0)

    # ---- Step 4: Marginal likelihood evaluation (C++) ----
    t0 <- proc.time()[3]
    ml_res <- marginal_likelihood_rb(
      XtX_cube = design$XtX,
      Xty_mat  = design$Xty,
      yty_vec  = design$yty,
      eta0     = eta0,
      V0       = V0,
      nu0      = nu0,
      psi0     = psi0,
      L        = L
    )
    log_ml <- as.numeric(ml_res$log_marginal_lik)
    timing$marglik <- timing$marglik + (proc.time()[3] - t0)

    # ---- Step 5: Weighting ----
    t0 <- proc.time()[3]
    log_target <- log_ml + log_prior_theta_const
    log_w <- log_target - draws$log_density
    log_w[!is.finite(log_w)] <- -Inf

    w <- exp(log_w - max(log_w))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)
    timing$weight <- timing$weight + (proc.time()[3] - t0)

    if (verbose >= 2) {
      cat(sprintf("\n\nRao-Blackwellized AIS - iter %d\n", iter))
      cat(sprintf("\n ESS = %.1f\n", ess))
    }

    if (ess >= min_ess) break
    if (iter == N_iter_max) {
      warning("maximum number of RB-AIS iterations reached")
      break
    }

    # ---- Step 5 (cont.): Update theta proposal ----
    t0 <- proc.time()[3]
    theta_prop_params <- update_theta_only_proposal(
      theta_unc = draws$theta_unc,
      w = w,
      prev_params = theta_prop_params,
      lr = lr
    )
    theta_prop_params$df <- theta_prop_params$df + eta_df_incr_per_iter
    if (iter >= first_iter_mult_N) {
      N_draw <- min(as.integer(N_draw * N_draw_mult), N_draw_max)
    }
    if (iter > 1 && ess < 0.8 * prev_ess) {
      theta_prop_params$Sigma <- theta_prop_params$Sigma * factor_inflate_Sigma
    }
    timing$update <- timing$update + (proc.time()[3] - t0)
    prev_ess <- ess
  }

  # ---- Step 6: Posterior Reconstruction ----
  t0 <- proc.time()[3]
  res_idx <- sample(N_draw, size = N_final, replace = TRUE, prob = w)
  thetas <- draws$theta[res_idx, , drop = FALSE]

  # For each resampled theta, draw sigma^2 and then eta
  posterior_scale <- as.numeric(ml_res$posterior_scale)
  nu_n <- nu0 + L

  sigma2s <- posterior_scale[res_idx] / stats::rchisq(N_final, df = nu_n)

  # Draw eta from MVN(mu_n, sigma^2 * Vn) using Cholesky
  etas <- matrix(0, nrow = N_final, ncol = n_eta)
  for (j in seq_len(N_final)) {
    idx <- res_idx[j]
    mu_n_j <- ml_res$mu_n[, idx]
    Vn_j <- ml_res$Vn[, , idx]
    # Cholesky of Vn (upper triangular), then L = t(R)
    R <- tryCatch(chol(Vn_j), error = function(e) {
      # If Cholesky fails, add small ridge and retry
      chol(Vn_j + 1e-8 * diag(n_eta))
    })
    z <- stats::rnorm(n_eta)
    etas[j, ] <- mu_n_j + sqrt(sigma2s[j]) * crossprod(R, z)
  }

  # Name the eta columns (C++ buffer order: s_m, s_{m-1}, ..., s_1)
  eta_col_names_cpp <- c("l", if (trend) "b",
                         if (seas) paste0("s", rev(seq_len(m))))
  colnames(etas) <- eta_col_names_cpp

  # Reorder to R convention (l, [b,] s1, ..., s_m) for RSS_vect_arma
  eta_col_names_r <- c("l", if (trend) "b",
                       if (seas) paste0("s", seq_len(m)))
  etas <- etas[, eta_col_names_r, drop = FALSE]

  # Compute final states by running the model forward with drawn (theta, eta)
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

  log_evidence <- max(log_w) + log(mean(exp(log_w - max(log_w))))

  log_lik_pointwise <- NULL
  if (return_pointwise) {
    E <- refit_final$residuals
    sd_mat <- matrix(sqrt(sigma2s), nrow = nrow(E), ncol = ncol(E), byrow = FALSE)
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
    prop_params = theta_prop_params,
    log_evidence = log_evidence,
    log_lik_pointwise = log_lik_pointwise,
    timing = timing
  )
}


##############################################################################
### Theta-Only Proposal Helpers ###

# Draw from a theta-only MVT proposal (no eta)
draw_theta_only <- function(N, proposal_params, theta_names,
                            phi_min, phi_max, antithetic = TRUE) {
  df <- proposal_params$df
  mu <- proposal_params$mus
  Sigma <- proposal_params$Sigma

  if (antithetic) {
    base <- mvtnorm::rmvt(n = ceiling(N / 2), sigma = Sigma, df = df, type = "shifted")
    samps_unc <- rbind(base, -base) +
      matrix(mu, nrow = nrow(base) * 2, ncol = length(mu), byrow = TRUE)
    samps_unc <- samps_unc[1:N, , drop = FALSE]
  } else {
    samps_unc <- mvtnorm::rmvt(n = N, sigma = Sigma, df = df,
                               delta = mu, type = "shifted")
  }
  colnames(samps_unc) <- theta_names

  log_density_unc <- mvtnorm::dmvt(samps_unc, delta = mu, sigma = Sigma,
                                   df = df, log = TRUE)

  # Transform to constrained space

  trans_res <- transform_unconstrained_to_theta(samps_unc, theta_names,
                                                phi_min, phi_max)
  log_density <- log_density_unc - trans_res$log_jac

  list(
    theta = trans_res$theta,
    theta_unc = samps_unc,
    log_density = log_density
  )
}

# Update theta-only proposal via weighted moment matching
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
