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
#' @param sampler Importance sampling algorithm: `"ais"` (single MVT proposal,
#'   default) or `"amis"` (Adaptive Multiple Importance Sampling with a
#'   mixture of MVT proposals and recycled historical samples).
#' @param rao_blackwellize_eta Logical; if `TRUE`, analytically integrate out
#'   the initial states (Rao-Blackwellization) instead of sampling them.
#'   Only supported with `sampler = "ais"`. Default is `FALSE`.
#' @param control Named list of tuning parameters. Missing values are filled from
#'   package defaults. NMIG-specific entries: `v_spike`, `v_slab`, `w_nmig`.
#'   `prior_models`: optional prior model probabilities for BMA.
#'   AMIS-specific entries: `K_mix` (number of mixture components, default 3),
#'   `ridge_eps` (ridge regularization, default 1e-4),
#'   `em_iter` (EM iterations per update, default 5),
#'   `jitter_scale` (initialization jitter, default 0.5).
#'
#' @return An object of class `"bets"`.
#' @export
bets <- function(y,
                 model = "ZZZ",
                 method = c("bma", "stacking", "nmig"),
                 additive.only = TRUE,
                 sampler = c("ais", "amis"),
                 rao_blackwellize_eta = FALSE,
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
  sampler <- match.arg(sampler)

  if (!is.logical(rao_blackwellize_eta) || length(rao_blackwellize_eta) != 1 ||
      is.na(rao_blackwellize_eta)) {
    stop("rao_blackwellize_eta must be TRUE or FALSE")
  }
  if (rao_blackwellize_eta && sampler == "amis") {
    stop("Rao-Blackwellized AMIS is not implemented. Please use AIS.")
  }

  ctrl <- resolve_bets_control(control)
  model_components <- coerce_model_components(model, stats::frequency(y), additive.only = TRUE)  # multiplicative models not supported yet

  fit <- fit_bets_models(
    y = y,
    model_components = model_components,
    ctrl = ctrl,
    method = method,
    sampler = sampler,
    rao_blackwellize_eta = rao_blackwellize_eta
  )

  structure(
    list(
      y = y,
      fit = fit,
      model_components = model_components,
      method = method,
      sampler = sampler,
      rao_blackwellize_eta = rao_blackwellize_eta,
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
  cat(sprintf("  sampler: %s\n", x$sampler))
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
