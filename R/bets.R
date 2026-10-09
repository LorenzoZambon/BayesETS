#' Fit a Bayesian ETS model
#'
#' Several ETS models are estimated separately in a Bayesian framework;
#' a posterior distribution is obtained for the parameters of each model.
#' The smoothing parameters are integrated by adaptive Gauss-Hermite quadrature (AGHQ)
#' or adaptive importance sampling (AIS);
#' the initial states and the error variance are integrated analytically.
#' Then the models are combined into a single predictive
#' distribution via Bayesian Model Averaging (BMA).
#' In this way, both parameter and model uncertainty are accounted for.
#' Currently only additive-error models are supported.
#'
#' @param y Univariate time series.
#' @param model Models to fit: either an ETS code (made of 3 or 4 characters, see **Details**)
#' or a vector of ETS codes. Default is `"ZZZ"` (all supported models).
#'   Duplicate models are dropped.
#' @param combination How the models are combined; only `"bma"` (default) is
#'   currently supported.
#' @param additive.only Logical; only `TRUE` (default) is currently supported.
#' @param verbose `0` (default) prints nothing,
#' `1` prints the model weights,
#' `2` also prints the progress of each fit.
#' `TRUE`/`FALSE` are accepted as `1`/`0`.
#' @param control Named list of tuning parameters, see **Details**. Missing
#'   entries take their default values.
#'
#' @details
#'
#' ## ETS codes
#'
#' An ETS code must have 3 or 4 characters; the possible values for each component are:
#' - Error: `"A"`, `"M"`, `"Z"`
#' - Trend: `"N"`, `"A"`, `"M"`, `"Z"`, `"Ad"`, `"Md"`
#' - Season: `"N"`, `"A"`, `"M"`, `"Z"`
#'
#' `"A"`: additive,
#' `"Ad"`: additive damped,
#' `"M"`: multiplicative,
#' `"Md"`: multiplicative damped,
#' `"Z"`: all supported options for that component.
#'
#' ## Integration methods
#'
#' Both methods start from the posterior mode of the smoothing
#' parameters (mapped to unconstrained space) and the inverse Hessian there.
#' Quadrature uses a Gauss-Hermite grid scaled by the inverse Hessian.
#' AIS samples a Student-t proposal and adapts it until the effective sample size (ESS)
#' reaches `min_ess`.
#' By default (`integration = "auto"`), quadrature is used
#' for models with up to 2 smoothing parameters (ANN, AAN, ANA) and AIS for
#' the others (AAdN, AAA, AAdA).
#'
#' ## Seasonal period
#'
#' The seasonal period is `frequency(y)`.
#' Seasonal models need an integer period between 2 and 24, and more than one full period of observations;
#' otherwise, only non-seasonal models can be fitted.
#'
#' ## Constant series
#'
#' A constant series is fitted with ETS(A,N,N), with a warning. The forecasts
#' are constant, and the width of the intervals depends on the prior of the
#' error variance (set by a heuristic, or by `psi0` in `control`).
#'
#' ## Control parameters
#'
#' The `control` argument accepts a named list with the following entries.
#' Entries indexed by `d` have one value per number of smoothing parameters
#' (`d` = 1, ..., 4); a scalar applies to all `d`.
#'
#' Integration:
#' - `integration`: `"auto"` (default), `"quadrature"` or `"ais"`.
#'   Quadrature is not recommended for models with 3 or
#'   more smoothing parameters.
#' - `n_scan`: number of prior points evaluated to start the search of the
#'   posterior mode. Default 64.
#' - `N_draw`: number of AIS draws per iteration, indexed by `d`. Default
#'   `c(128, 256, 512, 1024)`. Powers of 2 work best with Sobol points.
#' - `min_ess`, `N_iter_max`: AIS stops when the effective sample size reaches
#'   `min_ess` (default `N_draw / 4`), or after `N_iter_max` iterations
#'   (default 30). In the second case the model gets zero weight; if all
#'   models do, `bets()` stops with an error.
#' - `is_df`, `is_scale`: degrees of freedom (default 5) and scale inflation
#'   (default 4) of the AIS proposal.
#' - `lr`: maximum weight of the new estimates when the AIS proposal is
#'   updated (default 0.9).
#' - `n_quad`: number of Gauss-Hermite nodes per dimension, indexed by `d`.
#'   Default `c(21, 21, 9, 7)`.
#' - `N_final`: number of posterior draws kept per model, used by
#'   [predict.bets()], indexed by `d`. Default `c(200, 300, 500, 500)`.
#'
#' Priors:
#' - `nu0`, `psi0`: parameters of the prior of the error variance \eqn{\sigma^2},
#'   a scaled inverse chi-squared with `nu0` degrees of freedom (default 3, must be >2)
#'   and mean `psi0 / (nu0 - 2)`. By default, `psi0` is the
#'   residual variance of the naive or seasonal naive forecasts.
#' - `c_inflate_eta`: inflation factor of the prior covariance of the initial
#'   states, which is set by a heuristic (default 3).
#' - `phi_min`, `phi_max`: range of the damping parameter (default 0.8 and 0.98).
#'   All smoothing parameters have uniform priors.
#' - `prior_models`: prior model probabilities for BMA, one per model.
#'   Default: equal.
#'
#' @return An object of class `"bets"`: a list with the series `y`, the fit of
#'   each model and their weights (`fit`), the models (`model_components`), the
#'   `control` settings as given (with defaults), the value of `psi0` used, the
#'   seasonal `period` used for fitting, and the `call`.
#'
#' @examples
#' set.seed(1)
#' # Annual series: only additive non-seasonal models (ANN, AAN, AAdN), combined by BMA
#' fit <- bets(Nile, additive.only = TRUE)
#' fit
#'
#' # Monthly series: all additive models (ANN, AAN, AAdN, ANA, AAA, AAdA)
#' fit_m <- bets(USAccDeaths, additive.only = TRUE)
#' fit_m
#'
#' # Fit only some models
#' bets(USAccDeaths, model = c("ANA", "AAA"))
#'
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

  if (!is.numeric(y) || NCOL(y) != 1) {
    stop("y must be a numeric vector or a univariate time series")
  }
  if (any(!is.finite(y))) {
    stop("y must not contain missing or infinite values")
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
  if (combination == "stacking") {
    stop("Currently only Bayesian model averaging is supported (combination = \"bma\")")
  }

  if (!(is.numeric(verbose) || is.logical(verbose)) || length(verbose) != 1 || is.na(verbose)) {
    stop("verbose must be 0, 1, 2, TRUE or FALSE")
  }

  m <- seasonal_period(y)

  ctrl <- resolve_bets_control(control)
  model_components <- coerce_model_components(model, m, additive.only, n = length(y))
  if (!is.null(ctrl$prior_models) && length(ctrl$prior_models) != length(model_components)) {
    stop(sprintf("control$prior_models must have one value per model (%d)",
                 length(model_components)), call. = FALSE)
  }

  # Fit on a copy of y with frequency = seasonal period (1 if no seasonal model can be fitted);
  # y keeps its time index for predict()
  if (m > 24 || length(y) <= m) m <- 1L
  y_fit <- stats::ts(as.numeric(y), frequency = m)

  # Scale of the prior of sigma^2: data-based unless set in control (which keeps
  # the settings as given)
  psi0 <- if (is.null(ctrl$psi0)) default_psi0(y_fit) else ctrl$psi0
  fit_ctrl <- ctrl
  fit_ctrl$psi0 <- psi0

  if (isTRUE(all(y_fit == y_fit[1]))) {         # constant series
    warning("y is constant! Only ETS(A,N,N) is used: the forecasts are flat and ",
            "the prediction intervals might be unreliable.",
            call. = FALSE)
    fit <- fit_constant_series(y_fit, fit_ctrl, combination)
    model_components <- list(fit$results[[1]]$model_components)

  } else {                                      # non-constant series
    fit <- fit_bets_models(
      y = y_fit,
      model_components = model_components,
      ctrl = fit_ctrl,
      combination = combination,
      integration = ctrl$integration,
      verbose = as.numeric(verbose)
    )
  }

  structure(
    list(
      y = y,
      fit = fit,
      model_components = model_components,
      control = ctrl,
      psi0 = psi0,
      period = m,
      call = match.call()
    ),
    class = "bets"
  )
}

#' @export
print.bets <- function(x, ...) {
  cat("BETS model fit\n")
  cat(sprintf("  call: %s\n", paste(deparse(x$call, width.cutoff = 500L), collapse = " ")))
  cat(sprintf("  length of the series: %d\n", length(x$y)))
  cat(sprintf("  frequency: %s\n", format(stats::frequency(x$y))))

  print_comb <- ifelse(x$fit$combination == "bma", "Bayesian Model Averaging", "Stacking")
  cat(sprintf("  combination: %s\n", print_comb))

  # Models with weight >= 0.001, by decreasing weight; the others in one line
  mc <- x$model_components
  w <- x$fit$model_weights
  labels <- vapply(mc, ets_label, character(1))
  shown <- shown_models(w)
  cat("\n  Model weights and posterior means of the parameters:\n\n")
  cat_table(posterior_means_table(x$fit$results[shown], labels[shown], w[shown]))
  cat_negligible(labels, shown)

  if (x$fit$combination == "bma") cat_probabilities(component_probabilities(mc, w))
  invisible(x)
}

# Models (rows) with their weight and the posterior means of the smoothing
# parameters and of sigma ("-" if not in the model), as a character matrix
posterior_means_table <- function(results, labels, weights) {
  means <- lapply(results, function(r) {
    s <- posterior_summary(r$particles)
    stats::setNames(s$mean, rownames(s))
  })
  pars <- intersect(c("alpha", "beta", "phi", "gamma"), unlist(lapply(means, names)))
  par_cols <- lapply(pars, function(p) {
    vapply(means, function(m) if (p %in% names(m)) sprintf("%.3f", m[[p]]) else "-", character(1))
  })
  sigma <- vapply(means, `[[`, numeric(1), "sigma")
  tab <- cbind(labels, sprintf("%.3f", weights), do.call(cbind, par_cols), format(signif(sigma, 4)))
  colnames(tab) <- c("model", "weight", pars, "sigma")
  tab
}
