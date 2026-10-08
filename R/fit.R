################################################################################
# MODEL FITTING

# Fits each model and combines them (BMA or stacking)
fit_bets_models <- function(y,
                            model_components,
                            ctrl,
                            combination = c("bma", "stacking"),
                            integration = c("auto", "quadrature", "ais"),
                            verbose = 0) {
  combination <- match.arg(combination)
  integration <- match.arg(integration)

  ctrl$verbose <- verbose   # read by the integrators

  freq <- stats::frequency(y)

  n_models <- length(model_components)
  need_pointwise <- (combination == "stacking")

  # Centre the series at its initial level (only the level states change)
  y_shift <- mean(y[seq_len(max(1L, min(as.integer(freq), length(y))))])
  y_centred <- y - y_shift

  results_list <- vector("list", n_models)
  log_marginal_liks <- rep(NA_real_, n_models)
  log_lik_list <- if (need_pointwise) vector("list", n_models) else NULL
  fit_time_per_model <- numeric(n_models)

  for (i in seq_along(model_components)) {
    if (verbose >= 2) cat(sprintf("\nFitting model %d of %d\n", i, n_models))
    t0 <- proc.time()[3]

    integration_i <- resolve_integration(integration, model_components[[i]])
    res_i <- fit_one_model(y_centred, model_components[[i]], ctrl, integration_i,
                           return_pointwise = need_pointwise)
    # Undo the centring (NULL for failed models)
    if (!is.null(res_i$etas)) {
      res_i$etas[, "l"]   <- res_i$etas[, "l"] + y_shift
      res_i$states[, "l"] <- res_i$states[, "l"] + y_shift
    }

    fit_time_per_model[i] <- proc.time()[3] - t0
    results_list[[i]] <- res_i
    results_list[[i]]$model_components <- model_components[[i]]
    results_list[[i]]$integration <- integration_i
    log_marginal_liks[i] <- res_i$log_evidence
    if (need_pointwise) log_lik_list[[i]] <- res_i$log_lik_pointwise
  }

  if (!any(is.finite(log_marginal_liks))) {
    stop("No model could be fitted: all of them got zero weight (see the warnings). ",
         "Larger N_draw or N_iter_max in control may help AIS.", call. = FALSE)
  }

  t0 <- proc.time()[3]
  if (need_pointwise) {
    model_weights <- compute_stacking_weights(log_lik_list)
    if (verbose >= 1) cat("\nStacking Weights:\n")
  } else {
    prior_models <- ctrl$prior_models
    if (is.null(prior_models)) prior_models <- rep(1 / n_models, n_models)
    log_post_unnorm <- log(prior_models) + log_marginal_liks
    w <- exp(log_post_unnorm - max(log_post_unnorm))
    model_weights <- w / sum(w)
  }

  if (verbose >= 1) {
    cat("\nModel Weights:\n")
    labels <- vapply(model_components, ets_label, character(1))
    for (i in seq_along(labels)) {
      cat(sprintf("  %-5s: %.3f\n", labels[i], model_weights[i]))
    }
  }

  elapsed_combination <- proc.time()[3] - t0

  list(
    results = results_list,
    model_weights = model_weights,
    combination = combination,
    fit_time_per_model = fit_time_per_model,
    elapsed_combination = elapsed_combination
  )
}

# Fit of one model: prior, integrand, integration over the smoothing parameters
# (integrate_quadrature() or integrate_ais()), posterior draws
fit_one_model <- function(y, model_components, ctrl, integration,
                          return_pointwise = FALSE) {
  theta_names <- theta_names_of(model_components)
  N_final <- resolve_by_d(ctrl$N_final, length(theta_names))

  prior <- init_rb_prior(y, model_components, theta_names, ctrl)
  log_g_rb <- make_log_g_rb(y, model_components, theta_names, ctrl, prior)
  # Only log g, for the prior scan and the mode search
  log_g_value <- make_log_g_rb(y, model_components, theta_names, ctrl, prior,
                               posterior = FALSE)

  # The mode search starts at the best point of a prior scan
  t0 <- proc.time()[3]
  scan <- prior_scan(log_g_value, theta_names, ctrl$n_scan)
  t_scan <- proc.time()[3] - t0
  integrate <- switch(integration, quadrature = integrate_quadrature, ais = integrate_ais)
  int <- integrate(log_g_rb, scan$Z[1, ], theta_names, ctrl, log_g_mode = log_g_value)
  timing <- c(list(scan = t_scan), int$timing, list(post = 0))

  failed <- !is.null(int$failure)
  if (failed) {
    # The model gets zero weight; -1e300 rather than -Inf, to avoid NaN in logsumexp
    warning(int$failure, call. = FALSE)
    post <- list(log_lik_pointwise = if (return_pointwise) {
      matrix(-1e300, nrow = N_final, ncol = length(y))
    })
  } else {
    t0 <- proc.time()[3]
    post <- draw_rb_posterior(y, model_components, int$theta, int$w, int$ml_res,
                              N_final = N_final, nu0 = ctrl$nu0,
                              return_pointwise = return_pointwise)
    timing$post <- proc.time()[3] - t0
  }

  list(
    thetas = post$thetas,
    etas = post$etas,
    states = post$states,
    sigma2s = post$sigma2s,
    ess = int$ess,
    n_iter = int$n_iter,
    prop_params = int$proposal,
    log_evidence = if (failed) -Inf else int$log_evidence,
    log_lik_pointwise = post$log_lik_pointwise,
    timing = timing
  )
}

# Posterior draws: resampling of the weighted particles, then \sigma^2 and eta
# from their conditional posteriors, and the final states. ml_res: as returned
# by make_log_g_rb() at the particles.
draw_rb_posterior <- function(y, model_components, theta_particles, w, ml_res,
                              N_final, nu0, return_pointwise = FALSE) {
  L <- length(y)
  m <- stats::frequency(y)
  flags <- model_flags(model_components)
  n_eta <- n_states(model_components, m)
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
  eta_cpp <- ml_res$mu_n[, res_idx, drop = FALSE] +
    sweep(colSums(Rn_draws * Z_draws), 2, sqrt(sigma2s), "*")   # n_eta x N_final

  # Final states: affine functions of eta (C++)
  states <- t(final_states_rb(ml_res$final_coef, ml_res$final_const, eta_cpp,
                              as.integer(res_idx)))

  # Initial states from the C++ order (l, [b,] s_m, ..., s1) to the R order
  # (l, [b,] s1, ..., s_m); in the final states s1 is the next one used
  state_names <- c("l", if (flags$trend) "b", if (flags$seas) paste0("s", seq_len(m)))
  etas <- t(eta_cpp)
  colnames(etas) <- c("l", if (flags$trend) "b", if (flags$seas) paste0("s", rev(seq_len(m))))
  etas <- etas[, state_names, drop = FALSE]
  colnames(states) <- state_names

  log_lik_pointwise <- NULL
  if (return_pointwise) {
    # In-sample residuals, for stacking
    E <- ets_residuals(yR = as.numeric(y), trend = flags$trend, seas = flags$seas,
                       damped = flags$damped, m = m, init_statesR = etas,
                       paramsR = thetas, return_residuals = TRUE)$residuals
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

# Fit of a constant series with ETS(A,N,N).
# The posterior is simple: since residuals are all zero, alpha keeps its prior,
# the final level is the constant and \sigma^2 is driven by its prior only.
fit_constant_series <- function(y, ctrl, combination) {
  level <- as.numeric(y[1])
  psi0 <- ctrl$psi0
  n <- resolve_by_d(ctrl$N_final, 1)
  result <- list(
    thetas = matrix(stats::runif(n), ncol = 1, dimnames = list(NULL, "alpha")),
    etas = matrix(level, nrow = n, ncol = 1, dimnames = list(NULL, "l")),
    states = matrix(level, nrow = n, ncol = 1, dimnames = list(NULL, "l")),
    sigma2s = psi0 / stats::rchisq(n, df = ctrl$nu0 + length(y)),
    log_evidence = NA_real_,
    model_components = c("A", "N", "N", "FALSE"),
    integration = "constant series"
  )
  list(results = list(result), model_weights = 1, combination = combination,
       fit_time_per_model = 0, elapsed_combination = 0)
}

# Integration method of a model:
# with "auto", quadrature for up to 2 smoothing parameters, AIS otherwise
resolve_integration <- function(integration, model_components) {
  if (integration != "auto") return(integration)
  if (length(theta_names_of(model_components)) <= 2) "quadrature" else "ais"
}

# Stacking weights, from a list with one matrix of pointwise log-likelihoods
# (draws x time) per model
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
