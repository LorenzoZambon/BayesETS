################################################################################
# MODEL FITTING

# Negligible weights: the particles holding together at most .NEGLIGIBLE_PARTICLES
# of the weight of their model are dropped, and the models with weight below
# .NEGLIGIBLE_MODELS do not enter the forecasts
.NEGLIGIBLE_PARTICLES <- 1e-5
.NEGLIGIBLE_MODELS <- 1e-4

# Fits each model and combines them by BMA (stacking is not available: it needs
# pointwise predictive densities for compute_stacking_weights())
fit_bets_models <- function(y,
                            model_components,
                            ctrl,
                            integration = c("auto", "quadrature", "ais"),
                            verbose = 0) {
  integration <- match.arg(integration)

  ctrl$verbose <- verbose   # read by the integrators

  n_models <- length(model_components)

  # Centre the series at its initial level (only the level states change)
  y_shift <- centring_shift(y)
  y_centred <- y - y_shift

  results_list <- vector("list", n_models)
  log_marginal_liks <- rep(NA_real_, n_models)
  fit_time_per_model <- numeric(n_models)

  for (i in seq_along(model_components)) {
    if (verbose >= 2) {
      cat(sprintf("\nFitting model %-5s (%d of %d)\n", ets_label(model_components[[i]]), i, n_models))
    }
    t0 <- proc.time()[3]

    integration_i <- resolve_integration(integration, model_components[[i]])
    res_i <- fit_one_model(y_centred, model_components[[i]], ctrl, integration_i)
    # Undo the centring (no particles for failed models)
    if (!is.null(res_i$particles)) {
      res_i$particles$final_mean["l", ] <- res_i$particles$final_mean["l", ] + y_shift
    }

    fit_time_per_model[i] <- proc.time()[3] - t0
    results_list[[i]] <- res_i
    results_list[[i]]$model_components <- model_components[[i]]
    results_list[[i]]$integration <- integration_i
    log_marginal_liks[i] <- res_i$log_evidence
  }

  if (!any(is.finite(log_marginal_liks))) {
    stop("No model could be fitted: all of them got zero weight (see the warnings). ",
         "Larger N_draw or N_iter_max in control may help AIS.", call. = FALSE)
  }

  t0 <- proc.time()[3]
  prior_models <- ctrl$prior_models
  if (is.null(prior_models)) prior_models <- rep(1 / n_models, n_models)
  log_post_unnorm <- log(prior_models) + log_marginal_liks
  w <- exp(log_post_unnorm - max(log_post_unnorm))
  model_weights <- w / sum(w)

  if (verbose >= 1) {
    cat("\nModel Weights:\n")
    labels <- vapply(model_components, ets_label, character(1))
    for (i in seq_along(labels)) {
      cat(sprintf("  %-5s: %.3f\n", labels[i], model_weights[i]))
    }
  }

  # Models with negligible weight: final states dropped (no forecasts)
  for (i in which(model_weights < .NEGLIGIBLE_MODELS)) {
    if (!is.null(results_list[[i]]$particles)) {
      results_list[[i]]$particles[c("final_mean", "final_R")] <- NULL
    }
  }

  elapsed_combination <- proc.time()[3] - t0

  list(
    results = results_list,
    model_weights = model_weights,
    fit_time_per_model = fit_time_per_model,
    elapsed_combination = elapsed_combination
  )
}

# Fit of one model: prior, integrand, integration over the smoothing parameters
# (integrate_quadrature() or integrate_ais()), weighted particles of the posterior
fit_one_model <- function(y, model_components, ctrl, integration) {
  theta_names <- theta_names_of(model_components)

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

  # The model gets zero weight
  failed <- !is.null(int$failure)
  if (failed) warning(int$failure, call. = FALSE)

  # Weighted particles of the posterior (quadrature nodes or AIS draws), without the
  # negligible ones: theta, weights, \sigma^2 | theta ~ sigma2_scale / \chi^2_{sigma2_df},
  # and final states x_T | theta, \sigma^2 ~ N(final_mean, \sigma^2 t(final_R) final_R)
  particles <- if (!failed) {
    t0 <- proc.time()[3]
    keep <- heavy_set(int$w, .NEGLIGIBLE_PARTICLES)
    fs <- final_states_posterior_rb(int$ml_res$final_coef, int$ml_res$final_const,
                                    int$ml_res$mu_n, int$ml_res$Rn, keep)
    state_names <- state_names_of(model_components, stats::frequency(y))
    timing$post <- proc.time()[3] - t0
    list(theta = int$theta[keep, , drop = FALSE], w = int$w[keep] / sum(int$w[keep]),
         sigma2_scale = as.numeric(int$ml_res$posterior_scale)[keep],
         sigma2_df = ctrl$nu0 + length(y),
         final_mean = matrix(fs$mean, ncol = length(keep), dimnames = list(state_names, NULL)),
         final_R = array(fs$R, dim = dim(fs$R), dimnames = list(NULL, state_names, NULL)))
  }

  list(
    particles = particles,
    ess = int$ess,
    n_iter = int$n_iter,
    prop_params = int$proposal,
    log_evidence = if (failed) -Inf else int$log_evidence,
    timing = timing
  )
}

# Shift that centres the series at its initial level (mean of the first period)
centring_shift <- function(y) {
  mean(y[seq_len(max(1L, min(as.integer(stats::frequency(y)), length(y))))])
}

# Fit of a constant series with ETS(A,N,N).
# The posterior is simple: since residuals are all zero, alpha keeps its prior,
# the final level is the constant and \sigma^2 is driven by its prior only.
fit_constant_series <- function(y, ctrl) {
  level <- as.numeric(y[1])
  psi0 <- ctrl$psi0
  n <- 200L
  result <- list(
    # alpha keeps its uniform prior: n equally weighted, equally spaced quantiles
    particles = list(theta = matrix(stats::ppoints(n), ncol = 1, dimnames = list(NULL, "alpha")),
                     w = rep(1 / n, n), sigma2_scale = rep(psi0, n),
                     sigma2_df = ctrl$nu0 + length(y),
                     final_mean = matrix(level, nrow = 1, ncol = n, dimnames = list("l", NULL)),
                     final_R = array(0, dim = c(1, 1, n), dimnames = list(NULL, "l", NULL))),
    log_evidence = NA_real_,
    model_components = c("A", "N", "N", "FALSE"),
    integration = "constant series"
  )
  list(results = list(result), model_weights = 1, fit_time_per_model = 0,
       elapsed_combination = 0)
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
