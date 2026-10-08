# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("bets() rejects y with fewer than 3 observations", {
  expect_error(bets(c(1, 2)), "at least 3")
})

test_that("bets() rejects missing, infinite and non-numeric values", {
  expect_error(bets(ts(c(1, 2, NA, 4, 5))), "missing or infinite values")
  expect_error(bets(ts(c(1, 2, Inf, 4, 5))), "missing or infinite values")
  expect_error(bets(ts(letters[1:5])), "numeric vector or a univariate time series")
  expect_error(bets(ts(matrix(rnorm(20), ncol = 2))), "numeric vector or a univariate time series")
})

test_that("bets() rejects additive.only = FALSE", {
  expect_error(bets(ts(1:10), additive.only = FALSE), "additive\\.only = TRUE")
})

test_that("bets() rejects additive.only = NA", {
  expect_error(bets(ts(1:10), additive.only = NA), "additive\\.only must be TRUE or FALSE")
})

test_that("bets() rejects unknown control keys", {
  expect_error(bets(ts(1:10), control = list(foo = 1)), "Unknown control entry 'foo'")
})

test_that("bets() stops with a clear error when no model can be fitted", {
  set.seed(1)
  expect_error(
    suppressWarnings(bets(ts(rnorm(30)), model = list("ANN", "AAdN"),
                          control = list(integration = "ais", min_ess = 1e6, N_iter_max = 1))),
    "No model could be fitted")
})

# ---------------------------------------------------------------------------
# Return value structure
# ---------------------------------------------------------------------------

test_that("bets() coerces non-ts input to ts", {
  set.seed(1)
  fit <- bets(rnorm(20), model = "ANN")
  expect_true(stats::is.ts(fit$y))
})

test_that("bets() returns object of class 'bets'", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_s3_class(fit, "bets")
})

test_that("bets() result contains expected top-level elements", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_named(fit, c("y", "fit", "model_components", "control", "psi0", "call"),
               ignore.order = TRUE)
})

test_that("bets() model_weights sum to 1 under BMA with single model", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
})

test_that("bets() model_weights sum to 1 under BMA with multiple models", {
  # Single-model BMA is trivially 1; this multi-model case catches missing /sum(w)
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = list("ANN", "AAN"))
  expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
})

test_that("bets() rejects stacking (not implemented yet)", {
  expect_error(bets(ts(rnorm(30)), combination = "stacking"),
               "only Bayesian model averaging is supported")
})

test_that("bets() default integration ('auto') uses quadrature for d <= 2 and AIS for d >= 3", {
  set.seed(1)
  y <- ts(50 + cumsum(rnorm(32)) + rep(c(3, -1, -4, 2), 8), frequency = 4)
  fit <- bets(y, model = list("ANN", "AAN", "ANA", "AAdN", "AAA", "AAdA"))
  used <- vapply(fit$fit$results, `[[`, character(1), "integration")
  expect_equal(used, c("quadrature", "quadrature", "quadrature", "ais", "ais", "ais"))
})

test_that("control integration forces one method for all models", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = list("ANN", "AAdN"), control = list(integration = "ais"))
  expect_equal(vapply(fit$fit$results, `[[`, character(1), "integration"), c("ais", "ais"))
})

test_that("bets() rejects invalid combination, integration and verbose", {
  expect_error(bets(ts(rnorm(20)), model = "ANN", combination = "avg"), "should be one of")
  expect_error(bets(ts(rnorm(20)), model = "ANN", control = list(integration = "mcmc")),
               "integration must be one of")
  expect_error(bets(ts(rnorm(20)), model = "ANN", verbose = "yes"), "verbose must be")
})

test_that("bets() requires one prior probability per model", {
  expect_error(bets(ts(rnorm(30)), control = list(prior_models = c(0.5, 0.5))),
               "one value per model \\(3\\)")
  set.seed(1)
  fit <- bets(ts(rnorm(30)), model = list("ANN", "AAN"), control = list(prior_models = c(1, 0)))
  expect_equal(fit$fit$model_weights, c(1, 0))
})

test_that("bets() works for series with frequency < 1 (decennial uspop)", {
  set.seed(1)
  expect_warning(fit <- bets(uspop), "Frequency below 1")
  expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
  expect_true(all(is.finite(predict(fit, h = 3)$mean)))
  expect_output(print(fit), "frequency: 0.1")
})

test_that("non-integer frequency: warning, non-seasonal models, time index kept", {
  set.seed(1)
  y <- ts(10 + cumsum(rnorm(80)), start = c(2020, 1), frequency = 365.25 / 7)
  expect_warning(fit <- bets(y), "Non-integer seasonal period")
  expect_true(all(vapply(fit$model_components, `[[`, character(1), 3) == "N"))
  fc <- predict(fit, h = 3)
  expect_equal(stats::tsp(fc$mean)[1], stats::tsp(y)[2] + 7 / 365.25, tolerance = 1e-10)
  expect_output(print(fit), "frequency: 52.17857")
  expect_error(suppressWarnings(bets(y, model = "ANA")), "integer > 1")
})

test_that("frequency > 24: 'ZZZ' drops seasonal models, explicit ones are an error", {
  set.seed(1)
  y <- ts(10 + cumsum(rnorm(100)), frequency = 48)
  expect_warning(fit <- bets(y), "> 24")
  expect_true(all(vapply(fit$model_components, `[[`, character(1), 3) == "N"))
  expect_error(bets(y, model = "ANA"), "not supported for frequency\\(y\\) > 24")
})

test_that("seasonal series shorter than one period: non-seasonal models only", {
  set.seed(1)
  y <- ts(10 + cumsum(rnorm(10)), frequency = 12)
  fit <- bets(y)
  expect_length(fit$model_components, 3)
  expect_true(all(is.finite(predict(fit, h = 3)$mean)))
  expect_error(bets(y, model = "ANA"), "more than frequency\\(y\\) observations")
})

test_that("seasonal series with less than two periods: trend models get valid fits", {
  # The prior mean of b0 used to read observations m+1..2m, NA when L < 2m.
  set.seed(1)
  y <- ts(10 + cumsum(rnorm(20)), frequency = 12)
  fit <- suppressWarnings(bets(y))
  le <- vapply(fit$fit$results, `[[`, numeric(1), "log_evidence")
  expect_true(all(is.finite(le)))
})

test_that("constant series: flat forecasts, intervals from the sigma^2 prior", {
  width95 <- function(fc) as.numeric(fc$upper[, "95"] - fc$lower[, "95"])
  for (level in c(5, 0, -2e6)) {
    set.seed(1)
    expect_warning(fit <- bets(ts(rep(level, 20), frequency = 4)), "y is constant")
    expect_length(fit$model_components, 1)
    expect_output(print(fit), "ANN")
    fc <- predict(fit, h = 8)
    expect_lt(max(abs(fc$mean - level)), 0.01 * max(abs(level), 1))
    w <- width95(fc)
    expect_true(all(w > 0))
    expect_gt(w[8], w[1])   # wider with the horizon
  }
  # Narrower for longer series; the scale follows control$psi0
  set.seed(1)
  w_short <- width95(predict(suppressWarnings(bets(ts(rep(5, 5)))), h = 1, n_traj = 4000))
  set.seed(1)
  w_long <- width95(predict(suppressWarnings(bets(ts(rep(5, 50)))), h = 1, n_traj = 4000))
  expect_lt(w_long, w_short)
  w_psi <- vapply(c(1, 100), function(psi0) {
    set.seed(1)
    width95(predict(suppressWarnings(bets(ts(rep(5, 50)), control = list(psi0 = psi0))),
                    h = 1, n_traj = 4000))
  }, numeric(1))
  expect_equal(w_psi[2] / w_psi[1], 10, tolerance = 1e-6)
})

test_that("bets() fit contains one result per model", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = list("ANN", "AAN"))
  expect_length(fit$fit$results, 2)
})

test_that("print.bets() runs without error and prints header and call", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_output(print(fit), "BETS model fit")
  expect_output(print(fit), "call: bets\\(y = ts\\(rnorm\\(20\\)\\), model = \"ANN\"\\)")
  expect_false(any(grepl("Posterior probability", capture.output(print(fit)))))   # one model
})

test_that("print.bets() shows posterior means, negligible models and component probabilities", {
  set.seed(1)
  fit <- bets(ts(rnorm(30)), model = c("ANN", "AAN"), control = list(prior_models = c(1, 0)))
  out <- capture.output(print(fit))
  expect_true(any(grepl("^  model +weight +alpha +sigma$", out)))
  alpha <- sprintf("%.3f", mean(fit$fit$results[[1]]$thetas[, "alpha"]))
  expect_true(any(grepl(paste0("^  ANN +1\\.000 +", alpha, " "), out)))
  expect_true(any(grepl("(weight < 0.001: AAN)", out, fixed = TRUE)))
  expect_true(any(grepl("Posterior probability of: trend 0.00", out, fixed = TRUE)))
})

test_that("bets() stores the psi0 used, and keeps control as given", {
  set.seed(1)
  y <- ts(cumsum(rnorm(30)))
  fit <- bets(y, model = "ANN")
  expect_equal(fit$psi0, default_psi0(y))
  expect_null(fit$control$psi0)
  fit <- bets(y, model = "ANN", control = list(psi0 = 2))
  expect_equal(fit$psi0, 2)
  expect_equal(fit$control$psi0, 2)
  expect_equal(suppressWarnings(bets(ts(rep(5, 10))))$psi0, (0.2 * 5)^2)
})

