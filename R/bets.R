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
#'   `"stacking"`, or `"nmig"` (spike-and-slab model selection).
#' - `sampler`: Importance sampling algorithm: `"ais"` (single MVT proposal,
#'   default) or `"amis"` (Adaptive Multiple Importance Sampling with a
#'   mixture of MVT proposals and recycled historical samples).
#' - `rao_blackwellize_eta`: Logical; if `TRUE`, analytically integrates out
#'   the initial states (Rao-Blackwellization) instead of sampling them.
#'   Only supported with `sampler = "ais"`. Defaults to `TRUE` for AIS fits
#'   and `FALSE` for AMIS fits, where Rao-Blackwellization is not implemented.
#' - NMIG-specific: `v_spike`, `v_slab`, `w_nmig`.
#' - `prior_models`: optional prior model probabilities for BMA.
#' - AMIS-specific: `K_mix` (number of mixture components, default 3),
#'   `ridge_eps` (ridge regularization, default 1e-4),
#'   `em_iter` (EM iterations per update, default 5),
#'   `jitter_scale` (initialization jitter, default 0.5).
#' - `init`: Proposal initialization strategy: `"random_search"` (default),
#'   `"heuristic"`, or `"mle"`. When `"mle"`, Nelder-Mead optimization is run first
#'   to find an approximate MLE for theta (analytically integrating out eta at
#'   each evaluation). When `"random_search"`, a Sobol low-discrepancy sequence
#'   is evaluated over the unconstrained parameter space in a single vectorized
#'   batch call, and the best candidate is used as the proposal centre. This is
#'   typically faster than `"mle"` for low-dimensional models because all
#'   candidates are evaluated in one C++ call to `build_design_and_c_batch`.
#' - `n_sobol`: Number of Sobol candidates for `init = "random_search"`. Default
#'   `NULL` auto-selects `2^(d+3)` (16 / 32 / 64 / 128 for d = 1..4). Powers
#'   of 2 are optimal for Sobol sequences.
#' - `mle_tol`: Relative convergence tolerance for the Nelder-Mead optimizer
#'   used when `init = "mle"`. A high value (default `1e-3`) means we converge
#'   only roughly — sufficient to get a good starting region.
#' - `mle_maxit`: Maximum number of Nelder-Mead iterations (default `500`).
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
    match.arg(control$method, c("bma", "stacking", "nmig"))
  } else {
    "bma"
  }

  sampler <- if (!is.null(control$sampler)) {
    match.arg(control$sampler, c("ais", "amis"))
  } else {
    "ais"
  }

  rao_blackwellize_eta <- if (!is.null(control$rao_blackwellize_eta)) {
    rb_eta <- control$rao_blackwellize_eta
    if (!is.logical(rb_eta) || length(rb_eta) != 1 || is.na(rb_eta)) {
      stop("control$rao_blackwellize_eta must be TRUE or FALSE")
    }
    rb_eta
  } else {
    sampler != "amis"
  }

  if (rao_blackwellize_eta && sampler == "amis") {
    stop("Rao-Blackwellized AMIS is not implemented. Please use AIS.")
  }

  ctrl <- resolve_bets_control(control, rao_blackwellize_eta, freq)
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
