#' BETS package
#'
#' Bayesian ETS models.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom stats predict fitted residuals
#' @importFrom qrng sobol
#' @importFrom statmod gauss.quad
#' @useDynLib BETS, .registration = TRUE
"_PACKAGE"

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

# Default psi0: variance of the naive or (if smaller) of the seasonal naive residuals,
# floored to stay positive for deterministic series; heuristic for constant series.
default_psi0 <- function(y) {
  if (all(y == y[1])) return((0.2 * (if (y[1] != 0) abs(y[1]) else 1))^2)
  m <- stats::frequency(y)
  psi0 <- stats::var(diff(y))
  if (m > 1) psi0 <- min(psi0, stats::var(diff(y, lag = m)), na.rm = TRUE)
  max(psi0, 1e-8 * stats::var(y))
}

# Integration method of a model: 
# with "auto", quadrature for up to 2 smoothing parameters, AIS otherwise
resolve_integration <- function(integration, model_components) {
  if (integration != "auto") return(integration)
  if (length(theta_names_of(model_components)) <= 2) "quadrature" else "ais"
}
