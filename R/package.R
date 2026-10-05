#' BETS package
#'
#' Bayesian ETS models with adaptive importance sampling.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom stats predict
#' @importFrom qrng sobol
#' @importFrom statmod gauss.quad
#' @useDynLib BETS, .registration = TRUE
"_PACKAGE"

##############################################################################
### Model Fitting Wrapper ###

fit_bets_models <- function(y,
                            model_components,
                            ctrl,
                            method = c("bma", "stacking")) {
  method <- match.arg(method)
  integration_method <- match.arg(ctrl$integration_method, c("ais", "quadrature", "laplace_is"))
  integrate_model <- switch(integration_method,
                            ais        = adaptive_is_rb,
                            quadrature = quadrature_rb,
                            laplace_is = laplace_is_rb)

  verbose <- ctrl$verbose
  psi0    <- ctrl$psi0

  freq <- stats::frequency(y)

  n_models <- length(model_components)
  need_pointwise <- (method == "stacking")

  # set psi0 as either the MSE of the naive (if frequency = 1) or
  # the average of the MSEs of the naive and seasonal naive (if frequency > 1)
  if (is.null(psi0)) {
    mse_naive <- mean(diff(y, lag = 1)^2)
    psi0 <- if (freq > 1) 0.5 * (mean(diff(y, lag = freq)^2) + mse_naive) else mse_naive
    ctrl$psi0 <- psi0
  }

  results_list <- vector("list", n_models)
  log_marginal_liks <- rep(NA_real_, n_models)
  log_lik_list <- if (need_pointwise) vector("list", n_models) else NULL
  fit_time_per_model <- numeric(n_models)

  for (i in seq_along(model_components)) {
    if (verbose >= 2) cat(sprintf("\nFitting model %d of %d\n", i, n_models))
    t0 <- proc.time()[3]

    res_i <- integrate_model(
      y,
      model_components[[i]],
      ctrl = ctrl,
      return_pointwise = need_pointwise
    )

    fit_time_per_model[i] <- proc.time()[3] - t0
    results_list[[i]] <- res_i
    results_list[[i]]$model_components <- model_components[[i]]
    log_marginal_liks[i] <- res_i$log_evidence
    if (need_pointwise) log_lik_list[[i]] <- res_i$log_lik_pointwise
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
    fit_time_per_model = fit_time_per_model,
    elapsed_combination = elapsed_combination
  )
}
