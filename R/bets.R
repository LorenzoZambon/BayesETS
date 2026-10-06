#' Fit a Bayesian ETS model
#'
#' Estimates multiple ETS variants by drawing samples from the posterior distribution
#' of the parameters, using adaptive Gauss-Hermite quadrature or adaptive
#' importance sampling (AIS).
#' The resulting model fits are then combined into a single predictive distribution,
#' either by Bayesian Model Averaging (BMA) or stacking.
#' Currently only additive-error models are supported (i.e. `additive.only = TRUE`).
#'
#' @param y Univariate time series.
#' @param model Model space. Use `"ZZZ"` (default) to search a predefined set,
#'   or pass a specific model specification (e.g., "AAdN").
#' @param combination How the model fits are combined: `"bma"` (Bayesian Model
#'   Averaging, default) or `"stacking"`.
#' @param additive.only Logical; if `TRUE` (default), only additive-error models are considered.
#' Currently, setting `additive.only = FALSE` will result in an error.
#' @param verbose `0` (default) prints nothing, `1` prints the model weights,
#'   `2` also prints the progress of each model fit. `TRUE`/`FALSE` are
#'   accepted as `1`/`0`.
#' @param control Named list of tuning parameters. Missing values are filled from
#'   package defaults. See **Details**.
#'
#' @details
#' ## Integration methods
#'
#' The smoothing parameters of each model are integrated out either by adaptive
#' Gauss-Hermite quadrature or by adaptive importance sampling (AIS). Both start
#' from the posterior mode of the unconstrained smoothing parameters and the
#' inverse Hessian there. Quadrature integrates on a Gauss-Hermite grid scaled
#' by the inverse Hessian. AIS draws from a Student-t proposal centred at the
#' mode, sampled with randomised Sobol points, and adapts it only while the
#' effective sample size is below `min_ess` (usually a single step suffices).
#' Initial states and error variance are integrated out analytically in both
#' cases. By default (`integration = "auto"` in `control`), quadrature is used
#' for models with up to 2 smoothing parameters (ANN, AAN, ANA) and AIS for
#' models with 3 or 4 (AAdN, AAA, AAdA), where quadrature is less accurate.
#'
#' ## Control parameters
#'
#' The `control` argument accepts a named list with the following entries: (TODO)
#'
#' - `prior_models`: optional prior model probabilities for BMA, passed as a
#' list (...). Default: equal weights.
#' - `integration`: `"auto"` (default, see above), `"quadrature"` or `"ais"`,
#'   to force one method for all models. Forcing quadrature is not recommended
#'   for models with 3 or more smoothing parameters: its error grows with the
#'   dimension and can be large for strongly non-Gaussian posteriors.
#' - `N_draw`: Number of AIS draws per iteration, indexed by the number `d` of
#'   smoothing parameters (1 to 4). Default `c(128, 256, 512, 1024)`; a scalar
#'   applies to all `d`. Powers of 2 are optimal for Sobol sequences.
#' - `min_ess`: AIS stops as soon as the effective sample size of all draws so
#'   far reaches `min_ess`. Default `NULL` resolves to `N_draw / 4`. Models that
#'   do not reach it within `N_iter_max` iterations (default 30) get zero weight.
#' - `N_final`: Number of posterior draws kept per model (used by [predict.bets()]),
#'   indexed by `d` like `N_draw`. Default `c(100, 300, 500, 500)`.
#' - `is_df`, `is_scale`: Degrees of freedom of the Student-t proposal of AIS
#'   (default 5) and inflation of its initial scale matrix relative to the
#'   inverse Hessian (default 4).
#' - `n_scan`: Number of points of the prior scan (randomised Sobol points,
#'   evaluated in one batch) whose best point starts the search for the
#'   posterior mode, for both methods. Default 64.
#' - `n_quad`: Number of Gauss-Hermite nodes per dimension for quadrature,
#'   indexed by `d` like `N_draw`; the grid has `n_quad[d]^d` nodes. Default
#'   `c(21, 21, 9, 7)`.
#'
#' @return An object of class `"bets"`.
#' @export
bets <- function(y,
                 model = "ZZZ",
                 combination = c("bma", "stacking"),
                 additive.only = TRUE,
                 verbose = 0,
                 control = list()) {
  if (!stats::is.ts(y)) {
    y <- stats::ts(y)
  }
  if (length(y) < 3) {
    stop("y must contain at least 3 observations")
  }

  if (!is.logical(additive.only) || length(additive.only) != 1 || is.na(additive.only)) {
    stop("additive.only must be TRUE or FALSE")
  } else if (!additive.only) {
    stop("Currently only additive-error models are supported (additive.only = TRUE)")
  }

  combination <- match.arg(combination)
  if (!(is.numeric(verbose) || is.logical(verbose)) || length(verbose) != 1 || is.na(verbose)) {
    stop("verbose must be 0, 1, 2, TRUE or FALSE")
  }

  freq <- stats::frequency(y)

  ctrl <- resolve_bets_control(control, freq)
  integration <- match.arg(ctrl$integration, c("auto", "quadrature", "ais"))
  model_components <- coerce_model_components(model, freq, additive.only)

  fit <- fit_bets_models(
    y = y,
    model_components = model_components,
    ctrl = ctrl,
    combination = combination,
    integration = integration,
    verbose = as.numeric(verbose)
  )

  structure(
    list(
      y = y,
      fit = fit,
      model_components = model_components,
      control = ctrl,
      call = match.call()
    ),
    class = "bets"
  )
}

#' @export
print.bets <- function(x, ...) {
  cat("BETS model fit\n")
  cat(sprintf("  length of the series: %d\n", length(x$y)))
  cat(sprintf("  frequency: %d\n", stats::frequency(x$y)))

  print_comb <- ifelse(x$fit$combination == "bma", "Bayesian Model Averaging", "Stacking")
  cat(sprintf("  combination: %s\n", print_comb))

  mc <- x$model_components
  labels <- vapply(mc, ets_label, character(1))
  integration <- vapply(x$fit$results, `[[`, character(1), "integration")
  cat("\n  Models and weights:\n")
  for (i in seq_along(labels)) {
    cat(sprintf("  %-5s weight: %.3f  (%s)\n", labels[i], x$fit$model_weights[i], integration[i]))
  }
  invisible(x)
}
