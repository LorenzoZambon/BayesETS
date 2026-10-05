##############################################################################
### Rao-Blackwellized Adaptive Importance Sampling ###
###
### Samples only theta, analytically integrating out eta (initial states)
### and sigma^2 using conjugate Normal-Inverse-Chi-Squared priors.
###
### The first proposal is a multivariate Student-t centred at the mode of the
### integrand in the unconstrained coordinates z and scaled by its inverse
### Hessian (the Laplace fit shared with quadrature_rb()), sampled with
### randomised Sobol points.  It is adapted by weighted moment matching only
### while the ESS is below min_ess; the draws of all iterations are pooled.
##############################################################################

adaptive_is_rb <- function(y, model_components, ctrl,
                           return_pointwise = FALSE) {
  N_iter_max  <- ctrl$N_iter_max
  N_draw_raw  <- ctrl$N_draw    # may be a scalar or length-4 vector; resolved per d below
  N_final_raw <- ctrl$N_final   # same
  min_ess_raw <- ctrl$min_ess   # NULL (-> N_draw[d] / 4) or scalar or vector
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

  # Resolve dimension-dependent scalars now that d = n_theta is known.
  N_draw  <- resolve_by_d(N_draw_raw,  n_theta)
  N_final <- resolve_by_d(N_final_raw, n_theta)
  min_ess <- if (is.null(min_ess_raw)) N_draw / 4 else resolve_by_d(min_ess_raw, n_theta)

  # ---- Set up eta prior and the integrand ----
  prior <- init_rb_prior(y, model_components, theta_names, ctrl)
  log_g_rb <- make_log_g_rb(y, model_components, theta_names, ctrl, prior)

  # ---- Adaptive importance sampling ----
  ais <- adaptive_importance_sampling(log_g_rb, heuristic_z_start(theta_names),
                                      n_draw = N_draw, min_ess = min_ess,
                                      df = is_df, scale = is_scale,
                                      n_iter_max = N_iter_max, lr = lr,
                                      verbose = verbose)
  timing <- c(ais$timing, list(post = 0))

  prop_params <- ais$proposal
  names(prop_params$mus) <- theta_names
  rownames(prop_params$Sigma) <- colnames(prop_params$Sigma) <- theta_names

  # Convergence check: if ESS did not reach the target, skip posterior
  # reconstruction entirely — this model will receive zero weight in BMA/stacking.
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
      log_lik_pointwise = if (return_pointwise) matrix(-1e300, nrow = N_final, ncol = L) else NULL,  # use -1e300 instead of -Inf to avoid NaN in logsumexp
      timing = timing
    ))
  }

  # ---- Posterior Reconstruction (from the draws of all iterations) ----
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


##############################################################################
### Shared Rao-Blackwell Helpers (used by adaptive_is_rb and quadrature_rb) ###

# Conjugate prior for the initial states eta and sigma^2, plus the constant
# log-density of the uniform theta prior.  Returns eta0 and V0 already in the
# C++ ordering expected by build_design_and_c_batch / marginal_likelihood_rb.
init_rb_prior <- function(y, model_components, theta_names, ctrl) {
  psi0          <- ctrl$psi0
  phi_min       <- ctrl$phi_min
  phi_max       <- ctrl$phi_max
  c_inflate_eta <- ctrl$c_inflate_eta

  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  n_theta <- length(theta_names)

  # Use the existing heuristic for the initial states
  eta_init <- init_eta_params(y, model_components)
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

  # Expand to include the m-th seasonal slot, dropping the standard ETS
  # sum-to-zero constraint s_1 + ... + s_m = 0.  The design matrix X does NOT
  # enforce the constraint: the recursion is invariant under the joint shift
  # (l_0, s_1, ..., s_m) -> (l_0 + c, s_1 - c, ..., s_m - c), so X'X has rank
  # n_eta - 1.  We use a diagonal prior on all m slots; the prior regularizes
  # the unidentified shift direction without changing the predictive
  # distribution, and keeps V0 (and hence V0^{-1}) positive definite, which
  # is required by the Woodbury-based marginal_likelihood_rb kernel.
  # Encoding the constraint via off-diagonal coupling would make V0 exactly
  # rank-deficient and break inv_sympd(V0).
  if (seas && m > 1) {
    n_free <- length(eta0_free)
    n_full <- n_eta
    s_indices_free <- grep("^s\\d+$", names(eta0_free))
    Sigma_full <- matrix(0, n_full, n_full)
    Sigma_full[1:n_free, 1:n_free] <- Sigma_heuristic_free
    # Independent prior for the m-th seasonal slot, with the same marginal scale as the others
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

  # Uniform theta prior on the admissible region: its log-density is constant.
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

# Posterior reconstruction from weighted theta particles: resample N_final
# particles with probabilities w, then draw sigma^2 and eta from their
# conditional posteriors.  ml_res must come from marginal_likelihood_rb(...,
# return_posterior = TRUE) evaluated at theta_particles (one column/slice per row).
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

  # For each resampled theta, draw sigma^2 and then eta
  posterior_scale <- as.numeric(ml_res$posterior_scale)
  nu_n <- nu0 + L

  sigma2s <- posterior_scale[res_idx] / stats::rchisq(N_final, df = nu_n)

  # Draw eta from MVN(mu_n, sigma^2 * Vn), all N_final draws at once.
  # Rn = L_M^{-1} (lower triangular), pre-computed in C++. The formula
  # crossprod(Rn, z) = t(Rn) %*% z gives draws with covariance t(Rn)*Rn = Vn.
  # Column j of Z holds the standard normals of draw j (same RNG stream as one
  # rnorm(n_eta) per draw); Z_draws[k, i, j] = Z[k, j], so that
  # colSums(Rn_draws * Z_draws)[i, j] = (t(Rn_j) %*% Z[, j])[i].
  Z <- matrix(stats::rnorm(n_eta * N_final), nrow = n_eta, ncol = N_final)
  Rn_draws <- ml_res$Rn[, , res_idx, drop = FALSE]
  Z_draws  <- array(Z[, rep(seq_len(N_final), each = n_eta)], dim = c(n_eta, n_eta, N_final))
  etas <- t(ml_res$mu_n[, res_idx, drop = FALSE] +
              sweep(colSums(Rn_draws * Z_draws), 2, sqrt(sigma2s), "*"))

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


##############################################################################
### Adaptive Importance Sampling Helpers ###

# Adaptive importance sampling of g = exp(log g) over R^d.
#   log_g_fn, z_start : as in adaptive_gh_quadrature().
#   n_draw            : draws per iteration.
#   min_ess           : stop as soon as the ESS of the pooled draws reaches it.
# The first proposal is a multivariate Student-t (df) centred at the mode of
# log g with scale matrix scale * H^{-1} (H the Hessian of -log g there, see
# laplace_mode()).  While ESS < min_ess, the proposal is updated by weighted
# moment matching (update_theta_only_proposal()) and new draws are added.  The
# draws of all iterations are pooled and weighted against the mixture of all
# proposals used so far (deterministic-mixture weights), so the ESS can only
# grow.  Draws come from randomised Sobol points (RQMC).  Returns the log of
# the IS estimate of the integral, the normalised weights of the pooled draws,
# the draws, the output of log_g_fn per iteration and the ESS.
adaptive_importance_sampling <- function(log_g_fn, z_start, n_draw, min_ess,
                                         df = 5, scale = 1.5,
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
    # ---- Draw and evaluate log g, all draws of the iteration in one batch ----
    t0 <- proc.time()[3]
    Z_new <- draw_t_rqmc(n_draw, proposal$mus, proposal$Sigma, proposal$df)
    draw_eval <- log_g_fn(Z_new)
    timing$draws <- timing$draws + (proc.time()[3] - t0)
    proposals[[iter]]  <- c(proposal, list(n = nrow(Z_new)))
    draw_evals[[iter]] <- draw_eval
    Z <- rbind(Z, Z_new)
    log_g <- c(log_g, draw_eval$log_g)

    # ---- Weights of the pooled draws ----
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

    # ---- Update the proposal ----
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

# n RQMC draws from the multivariate Student-t (df) with location mu and scale
# matrix Sigma: d coordinates of randomised Sobol points through qnorm, one
# through qchisq for the mixing variable.
draw_t_rqmc <- function(n, mu, Sigma, df) {
  d <- length(mu)
  u <- matrix(qrng::sobol(n, d = d + 1, randomize = "digital.shift"), ncol = d + 1)
  eps  <- stats::qnorm(u[, seq_len(d), drop = FALSE])
  chi2 <- stats::qchisq(u[, d + 1], df = df)
  devs <- (eps %*% chol(Sigma)) / sqrt(chi2 / df)
  # Drop the (measure-zero) draws from points that land exactly on 0 or 1.
  devs <- devs[rowSums(!is.finite(devs)) == 0, , drop = FALSE]
  sweep(devs, 2, mu, "+")
}

# Log density at the rows of Z of the mixture of Student-t proposals, each
# weighted by its share of draws (list elements: mus, Sigma, df, n).
log_mix_density <- function(Z, proposals) {
  n_tot <- sum(vapply(proposals, `[[`, numeric(1), "n"))
  log_q <- vapply(proposals, function(p) {
    log(p$n / n_tot) + ldmvt_chol(sweep(Z, 2, p$mus, "-"), chol(p$Sigma), p$df)
  }, numeric(nrow(Z)))
  log_q <- matrix(log_q, nrow = nrow(Z))
  q_max <- apply(log_q, 1, max)
  q_max + log(rowSums(exp(log_q - q_max)))
}

# Concatenate marginal_likelihood_rb(..., return_posterior = TRUE) outputs of
# several batches of particles, in the same order as their rows.
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
