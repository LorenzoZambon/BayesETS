##############################################################################
### Rao-Blackwellized Importance Sampling with a Laplace Proposal ###
###
### Alternative to adaptive_is_rb() and quadrature_rb(): a single importance
### sampling step in the unconstrained coordinates z, with a multivariate
### Student-t proposal centred at the mode of the integrand and scaled by its
### inverse Hessian (the same Laplace fit as quadrature_rb()), sampled with
### randomised Sobol points.  No adaptation.  eta and sigma^2 are integrated
### out analytically exactly as in adaptive_is_rb().
##############################################################################

laplace_is_rb <- function(y, model_components, ctrl,
                          return_pointwise = FALSE) {
  N_final_raw <- ctrl$N_final
  n_is_raw    <- ctrl$n_is      # scalar or length-4 vector; resolved per d below
  is_df       <- ctrl$is_df
  is_scale    <- ctrl$is_scale
  nu0         <- ctrl$nu0
  verbose     <- ctrl$verbose

  L <- length(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha", if (trend) c("beta", if (damped) "phi"), if (seas) "gamma")
  n_theta <- length(theta_names)

  # Resolve dimension-dependent scalars now that d = n_theta is known.
  N_final <- resolve_by_d(N_final_raw, n_theta)
  n_is    <- resolve_by_d(n_is_raw,    n_theta)

  # ---- Set up eta prior and the integrand ----
  prior <- init_rb_prior(y, model_components, theta_names, ctrl)
  log_g_rb <- make_log_g_rb(y, model_components, theta_names, ctrl, prior)

  # ---- Importance sampling ----
  lis <- laplace_importance_sampling(log_g_rb, heuristic_z_start(theta_names), n_is,
                                     df = is_df, scale = is_scale)
  timing <- c(lis$timing, list(post = 0))

  Sigma <- lis$Sigma
  rownames(Sigma) <- colnames(Sigma) <- theta_names
  prop_params <- list(mus = stats::setNames(lis$zhat, theta_names), Sigma = Sigma, df = is_df)

  if (verbose >= 2) {
    cat(sprintf("\n\nRB Laplace importance sampling: %d points\n", n_is))
    cat(sprintf("\n ESS = %.1f, log evidence = %.4f\n", lis$ess, lis$log_evidence))
  }

  # No finite weight: this model will receive zero weight in BMA/stacking.
  if (!is.finite(lis$log_evidence)) {
    warning("laplace_is_rb: log g is not finite at any draw")
    return(list(
      thetas = NULL,
      etas = NULL,
      states = NULL,
      sigma2s = NULL,
      ess = 0,
      n_iter = NA_integer_,
      prop_params = prop_params,
      log_evidence = -Inf,
      log_lik_pointwise = if (return_pointwise) matrix(-1e300, nrow = N_final, ncol = L) else NULL,  # use -1e300 instead of -Inf to avoid NaN in logsumexp
      timing = timing
    ))
  }
  if (lis$ess < 0.1 * n_is) {
    warning(sprintf(paste0(
      "laplace_is_rb: low effective sample size (%.0f of %d draws); ",
      "the Laplace proposal fits the posterior poorly"), lis$ess, n_is))
  }

  # ---- Posterior Reconstruction ----
  t0 <- proc.time()[3]
  post <- draw_rb_posterior(y, model_components,
                            theta_particles = lis$draw_eval$theta,
                            w = lis$w,
                            ml_res = lis$draw_eval$ml_res,
                            N_final = N_final, nu0 = nu0,
                            return_pointwise = return_pointwise)
  timing$post <- proc.time()[3] - t0

  list(
    thetas = post$thetas,
    etas = post$etas,
    states = post$states,
    sigma2s = post$sigma2s,
    ess = lis$ess,
    n_iter = NA_integer_,
    prop_params = prop_params,
    log_evidence = lis$log_evidence,
    log_lik_pointwise = post$log_lik_pointwise,
    timing = timing
  )
}


##############################################################################
### Laplace Importance Sampling Helper ###

# Importance sampling of g = exp(log g) over R^d with a multivariate Student-t
# proposal (df degrees of freedom) centred at the mode of log g, with scale
# matrix scale * H^{-1} (H the Hessian of -log g there, see laplace_mode()).
#   log_g_fn, z_start : as in adaptive_gh_quadrature().
#   n_is              : number of draws.
# The draws come from randomised Sobol points (RQMC): d coordinates through
# qnorm and one through qchisq for the mixing variable of the t.  Returns the
# log of the IS estimate of the integral, the normalised weights, the draws,
# the output of log_g_fn at the draws and the effective sample size.
laplace_importance_sampling <- function(log_g_fn, z_start, n_is,
                                        df = 5, scale = 1.5, ...) {
  d <- length(z_start)
  lap <- laplace_mode(log_g_fn, z_start, ...)
  Sigma <- scale * tcrossprod(lap$Lmat)
  R_chol <- chol(Sigma)
  timing <- c(lap$timing, list(draws = 0))

  # ---- RQMC draws from the t proposal ----
  u <- matrix(qrng::sobol(n_is, d = d + 1, randomize = "digital.shift"), ncol = d + 1)
  eps  <- stats::qnorm(u[, seq_len(d), drop = FALSE])
  chi2 <- stats::qchisq(u[, d + 1], df = df)
  devs <- (eps %*% R_chol) / sqrt(chi2 / df)
  # Drop the (measure-zero) draws from points that land exactly on 0 or 1.
  devs <- devs[rowSums(!is.finite(devs)) == 0, , drop = FALSE]
  Z <- sweep(devs, 2, lap$zhat, "+")

  # ---- Weights: g / q, all draws in one batch ----
  t0 <- proc.time()[3]
  draw_eval <- log_g_fn(Z)
  log_w <- draw_eval$log_g - ldmvt_chol(devs, R_chol, df)
  log_w[!is.finite(log_w)] <- -Inf
  timing$draws <- proc.time()[3] - t0

  lw_max <- max(log_w)
  if (is.finite(lw_max)) {
    w <- exp(log_w - lw_max)
    log_evidence <- lw_max + log(mean(w))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)
  } else {
    w <- NULL
    log_evidence <- -Inf
    ess <- 0
  }

  list(
    log_evidence = log_evidence,
    w = w,
    Z = Z,
    draw_eval = draw_eval,
    ess = ess,
    zhat = lap$zhat,
    Sigma = Sigma,
    optim = lap$optim,
    timing = timing
  )
}
