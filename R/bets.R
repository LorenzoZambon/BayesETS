#' Fit a Bayesian ETS model
#'
#' Fits additive-error ETS variants using adaptive importance sampling,
#' then combines models with BMA or stacking.
#'
#' @param y Univariate time series.
#' @param model Model space. Use `"ZZZ"` (default) to search a predefined set,
#'   or pass a specific model specification.
#' @param method Combination strategy for multi-model fitting:
#'   `"bma"` (Bayesian Model Average, default), `"stacking"`, or `"nmig"` (spike-and-slab model selection).
#' @param additive.only Logical; if `TRUE` (default), multiplicative ETS components are rejected.
#' @param control Named list of tuning parameters. Missing values are filled from
#'   package defaults. NMIG-specific entries: `v_spike`, `v_slab`, `w_nmig`.
#'   `prior_models`: optional prior model probabilities for BMA.
#'
#' @return An object of class `"bets"`.
#' @export
bets <- function(y,
                 model = "ZZZ",
                 method = c("bma", "stacking", "nmig"),
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

  method <- match.arg(method)
  ctrl <- resolve_bets_control(control)
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
      method = method,
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
  cat(sprintf("  combination: %s\n", x$method))

  mc <- if (x$method == "nmig") {
    lapply(x$fit$results, `[[`, "model_components")
  } else {
    x$model_components
  }
  labels <- vapply(mc, ets_label, character(1))
  for (i in seq_along(labels)) {
    cat(sprintf("  %-5s weight: %.3f\n", labels[i], x$fit$model_weights[i]))
  }
  invisible(x)
}
