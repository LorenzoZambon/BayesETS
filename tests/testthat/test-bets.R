# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("bets() rejects y with fewer than 3 observations", {
  expect_error(bets(c(1, 2)), "at least 3")
})

test_that("bets() rejects additive.only = FALSE", {
  expect_error(bets(ts(1:10), additive.only = FALSE), "additive\\.only = TRUE")
})

test_that("bets() rejects additive.only = NA", {
  expect_error(bets(ts(1:10), additive.only = NA), "additive\\.only must be TRUE or FALSE")
})

test_that("bets() rejects unknown control keys", {
  expect_error(bets(ts(1:10), control = list(foo = 1)), "Unknown control entries")
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
  expect_named(fit, c("y", "fit", "model_components", "control", "call"),
               ignore.order = TRUE)
})

test_that("bets() model_weights sum to 1 under BMA", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
})

test_that("bets() model_weights sum to 1 under stacking", {
  set.seed(1)
  fit <- bets(ts(rnorm(30)), model = list("ANN", "AAN"),
              control = list(method = "stacking"))
  expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
})

test_that("bets() fit contains one result per model", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = list("ANN", "AAN"))
  expect_length(fit$fit$results, 2)
})

test_that("print.bets() runs without error and prints header", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN")
  expect_output(print(fit), "BETS model fit")
})

