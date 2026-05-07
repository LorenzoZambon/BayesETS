#' Fit a Bayesian ETS model
#'
#' Estimates multiple ETS variants by drawing samples from the posterior distribution
#' of the parameters, using adaptive importance sampling (AIS).
#' The resulting model fits are then combined into a single predictive distribution,
#' either by Bayesian Model Averaging (BMA) or stacking.
#' Currently only additive-error models are supported (i.e. `additive.only = TRUE`).
#'
#' @param y Univariate time series.
#' @param model Model space. Use `"ZZZ"` (default) to search a predefined set,
#'   or pass a specific model specification (e.g., "AAdN").
#' @param additive.only Logical; if `TRUE` (default), only additive-error models are considered.
#' Currently, setting `additive.only = FALSE` will result in an error.
#' @param control Named list of tuning parameters. Missing values are filled from
#'   package defaults. See **Details**.
#'
#' @details
#' ## Control parameters
#'
#' The `control` argument accepts a named list with the following entries: (TODO)
#'
#' - `method`: Combination strategy: `"bma"` (Bayesian Model Average, default),
#'   or `"stacking"`.
#' - `prior_models`: optional prior model probabilities for BMA, passed as a
#' list (...). Default: equal weights.
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
    control$method <- NULL  # remove method from control to avoid confusion later
  } else {
    "bma"
  }

  ctrl <- resolve_bets_control(control, freq)
  model_components <- coerce_model_components(model, freq, additive.only)

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
  cat(sprintf("  length of the series: %d\n", length(x$y)))
  cat(sprintf("  frequency: %d\n", stats::frequency(x$y)))

  print_comb <- ifelse(x$control$method == "bma", "Bayesian Model Averaging", "Stacking")
  cat(sprintf("  combination: %s\n", print_comb))

  mc <- x$model_components
  labels <- vapply(mc, ets_label, character(1))
  cat("\n  Models and weights:\n")
  for (i in seq_along(labels)) {
    cat(sprintf("  %-5s weight: %.3f\n", labels[i], x$fit$model_weights[i]))
  }
  invisible(x)
}
