#' Fit a Bayesian ETS model
#'
#' Fits additive-error ETS variants using adaptive importance sampling,
#' then combines models with BMA or stacking.
#'
#' @param y Univariate time series.
#' @param model Model space. Use `"ZZZ"` (default) to search a predefined set,
#'   or pass a specific model specification.
#' @param additive.only Logical; if `TRUE` (default), multiplicative ETS components are rejected.
#' @param control Named list of tuning parameters. Missing values are filled from
#'   package defaults. See **Details**.
#'
#' @details
#' ## Control parameters
#'
#' The `control` argument accepts a named list with the following entries:
#'
#' - `method`: Combination strategy: `"bma"` (Bayesian Model Average, default),
#'   or `"stacking"`.
#' - `prior_models`: optional prior model probabilities for BMA.
#' - `n_sobol`: Number of Sobol candidates used to build the proposal distribution. 
#'   Default `NULL` resolves to `N_draw` (same budget
#'   as one AIS iteration), which already scales with problem dimension through
#'   `N_draw`. Powers of 2 are optimal for Sobol sequences.
#'
#' @return An object of class `"bets"`.
#' @export
bets <- function(y,
                 model = "ZZZ",
                 additive.only = TRUE,
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

  freq <- stats::frequency(y)

  method <- if (!is.null(control$method)) {
    match.arg(control$method, c("bma", "stacking"))
  } else {
    "bma"
  }

  ctrl <- resolve_bets_control(control, freq)
  model_components <- coerce_model_components(model, stats::frequency(y), additive.only = TRUE)  # multiplicative models not supported yet

  fit <- fit_bets_models(
    y = y,
    model_components = model_components,
    ctrl = ctrl,
    method = method
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
  cat(sprintf("  length(y): %d\n", length(x$y)))
  cat(sprintf("  frequency: %d\n", stats::frequency(x$y)))
  cat(sprintf("  combination: %s\n", x$control$method))

  mc <- x$model_components
  labels <- vapply(mc, ets_label, character(1))
  for (i in seq_along(labels)) {
    cat(sprintf("  %-5s weight: %.3f\n", labels[i], x$fit$model_weights[i]))
  }
  invisible(x)
}
