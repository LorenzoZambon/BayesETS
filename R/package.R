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
                            combination = c("bma", "stacking"),
                            integration = c("auto", "quadrature", "ais"),
                            verbose = 0) {
  combination <- match.arg(combination)
  integration <- match.arg(integration)

  ctrl$verbose <- verbose   # read by the integrators
  psi0    <- ctrl$psi0

  freq <- stats::frequency(y)

  n_models <- length(model_components)
  need_pointwise <- (combination == "stacking")

  # set psi0 as either the MSE of the naive (if frequency = 1) or
  # the average of the MSEs of the naive and seasonal naive (if frequency > 1)
  if (is.null(psi0)) {
    mse_naive <- mean(diff(y, lag = 1)^2)
    psi0 <- if (freq > 1) 0.5 * (mean(diff(y, lag = freq)^2) + mse_naive) else mse_naive
    ctrl$psi0 <- psi0
  }

  # Centre the series at its initial level (the prior mean of l0, see
  # init_eta_params(): mean of the first frequency observations, at least one,
  # since the frequency can be < 1, e.g. decennial data).  The additive model
  # is location-equivariant: only the level states shift, so the evidence and
  # theta posterior are unchanged, while the sufficient statistics of the C++
  # kernels stay small.
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
    integrate_model <- switch(integration_i,
                              ais        = adaptive_is_rb,
                              quadrature = quadrature_rb)
    res_i <- integrate_model(
      y_centred,
      model_components[[i]],
      ctrl = ctrl,
      return_pointwise = need_pointwise
    )
    # Back to the original location (NULL if the model failed)
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

# Integration method for one model: "auto" uses quadrature for up to 2
# smoothing parameters (most accurate and cheapest there) and AIS for 3-4
# (quadrature has a downward bias there that grows with the dimension).
resolve_integration <- function(integration, model_components) {
  if (integration != "auto") return(integration)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  n_theta <- 1 + (if (trend) 1 + damped else 0) + seas
  if (n_theta <= 2) "quadrature" else "ais"
}
