#' BETS package
#'
#' Bayesian ETS models.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom stats predict
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
  psi0    <- ctrl$psi0

  freq <- stats::frequency(y)

  n_models <- length(model_components)
  need_pointwise <- (combination == "stacking")

  if (is.null(psi0)) {
    psi0 <- default_psi0(y)
    ctrl$psi0 <- psi0
  }

  # Centre the series at its initial level: only the level states change,
  # and the C++ kernels work with smaller numbers
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

# Fit of a constant series with ETS(A,N,N). The integrators cannot be used (all
# residuals are zero), but the posterior is simple: alpha keeps its prior, the
# final level is the constant and \sigma^2 is driven by its prior only.
fit_constant_series <- function(y, ctrl, combination) {
  level <- as.numeric(y[1])
  psi0 <- ctrl$psi0
  if (is.null(psi0)) psi0 <- (0.2 * (if (level != 0) abs(level) else 1))^2
  n <- resolve_by_d(ctrl$N_final, 1)
  result <- list(
    thetas = matrix(stats::runif(n), ncol = 1, dimnames = list(NULL, "alpha")),
    states = matrix(level, nrow = n, ncol = 1, dimnames = list(NULL, "l")),
    sigma2s = psi0 / stats::rchisq(n, df = ctrl$nu0 + length(y)),
    log_evidence = NA_real_,
    model_components = c("A", "N", "N", "FALSE"),
    integration = "constant series"
  )
  list(results = list(result), model_weights = 1, combination = combination,
       fit_time_per_model = 0, elapsed_combination = 0)
}

# Default psi0: variance of the naive or, if smaller, of the seasonal naive
# residuals (variance rather than MSE, to ignore drift). Floored to stay
# positive for deterministic series.
default_psi0 <- function(y) {
  m <- stats::frequency(y)
  psi0 <- stats::var(diff(y))
  if (m > 1) psi0 <- min(psi0, stats::var(diff(y, lag = m)), na.rm = TRUE)
  max(psi0, 1e-8 * stats::var(y))
}

# Integration method of a model: with "auto", quadrature for up to 2 smoothing
# parameters, AIS otherwise
resolve_integration <- function(integration, model_components) {
  if (integration != "auto") return(integration)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  n_theta <- 1 + (if (trend) 1 + damped else 0) + seas
  if (n_theta <= 2) "quadrature" else "ais"
}
