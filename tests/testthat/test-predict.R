# Shared fitted objects, created once for the whole file
local({
  set.seed(42)
  fit_ann  <<- bets(ts(rnorm(20)), model = "ANN")

  set.seed(42)
  fit_seas <<- bets(ts(rnorm(24), frequency = 4), model = "ANA")
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("predict.bets() errors when newdata is supplied", {
  expect_error(predict(fit_ann, newdata = 1:5), "newdata")
})

test_that("predict.bets() errors for non-positive or non-integer h", {
  expect_error(predict(fit_ann, h = 0),   "h must be a single positive integer")
  expect_error(predict(fit_ann, h = -1),  "h must be a single positive integer")
  expect_error(predict(fit_ann, h = 0.5), "h must be a single positive integer")
  expect_error(predict(fit_ann, h = 2.7), "h must be a single positive integer")
  expect_error(predict(fit_ann, h = NA),  "h must be a single positive integer")
})

test_that("predict.bets() errors for non-scalar h", {
  expect_error(predict(fit_ann, h = c(5, 10)), "h must be a single positive integer")
})

test_that("predict.bets() errors for invalid n_traj", {
  expect_error(predict(fit_ann, h = 5, n_traj = 0), "n_traj must be a single positive integer")
  expect_error(predict(fit_ann, h = 5, n_traj = 10.5), "n_traj must be a single positive integer")
  expect_error(predict(fit_ann, h = 5, n_traj = c(100, 200)), "n_traj must be a single positive integer")
})

test_that("predict.bets() accepts a custom n_traj", {
  fc <- predict(fit_ann, h = 5, n_traj = 50)
  expect_true(all(is.finite(fc$mean)))
})

test_that("predict.bets() errors when level is out of (0, 100)", {
  expect_error(predict(fit_ann, h = 5, level = 0),   "level must be between 0 and 100")
  expect_error(predict(fit_ann, h = 5, level = 100), "level must be between 0 and 100")
})

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("predict.bets() returns class 'bets_forecast'", {
  fc <- predict(fit_ann, h = 5)
  expect_s3_class(fc, "bets_forecast")
})

test_that("predict.bets() returns expected list elements", {
  fc <- predict(fit_ann, h = 5)
  expect_named(fc, c("method", "model", "mean", "lower", "upper", "level", "x", "series"),
               ignore.order = TRUE)
})

# ---------------------------------------------------------------------------
# mean forecast dimensions and timing
# ---------------------------------------------------------------------------

test_that("predict.bets() mean is a ts of length h", {
  fc <- predict(fit_ann, h = 7)
  expect_true(stats::is.ts(fc$mean))
  expect_equal(length(fc$mean), 7)
})

test_that("predict.bets() mean starts immediately after end of training series", {
  fc <- predict(fit_ann, h = 5)
  expected_start <- stats::tsp(fit_ann$y)[2] + 1 / stats::frequency(fit_ann$y)
  expect_equal(stats::tsp(fc$mean)[1], expected_start, tolerance = 1e-10)
})

test_that("predict.bets() mean starts correctly for seasonal series (freq > 1)", {
  fc <- predict(fit_seas, h = 4)
  # Time step 1 / frequency
  expected_start <- stats::tsp(fit_seas$y)[2] + 1 / stats::frequency(fit_seas$y)
  expect_equal(stats::tsp(fc$mean)[1], expected_start, tolerance = 1e-10)
})

test_that("predict.bets() mean frequency matches training series", {
  fc_ann  <- predict(fit_ann,  h = 4)
  fc_seas <- predict(fit_seas, h = 4)
  expect_equal(stats::frequency(fc_ann$mean),  stats::frequency(fit_ann$y))
  expect_equal(stats::frequency(fc_seas$mean), stats::frequency(fit_seas$y))
})

# ---------------------------------------------------------------------------
# Interval dimensions
# ---------------------------------------------------------------------------

test_that("predict.bets() lower and upper have h rows", {
  fc <- predict(fit_ann, h = 6)
  expect_equal(nrow(fc$lower), 6)
  expect_equal(nrow(fc$upper), 6)
})

test_that("predict.bets() works for h = 1 and for a single level", {
  fc <- predict(fit_ann, h = 1)
  expect_equal(dim(fc$lower), c(1L, 2L))
  expect_equal(colnames(fc$upper), c("80", "95"))
  fc <- predict(fit_ann, h = 4, level = 90)
  expect_equal(dim(fc$lower), c(4L, 1L))
  fc <- predict(fit_ann, h = 1, level = 90)
  expect_equal(dim(fc$upper), c(1L, 1L))
})

test_that("predict.bets() lower and upper have ncol = length(level)", {
  fc <- predict(fit_ann, h = 5, level = c(80, 90, 95))
  expect_equal(ncol(fc$lower), 3)
  expect_equal(ncol(fc$upper), 3)
})

test_that("predict.bets() lower < upper at every step and level", {
  fc <- predict(fit_ann, h = 5)
  expect_true(all(fc$lower < fc$upper))
})

test_that("predict.bets() lower and upper column names match level", {
  fc <- predict(fit_ann, h = 5, level = c(80, 95))
  expect_equal(colnames(fc$lower), c("80", "95"))
  expect_equal(colnames(fc$upper), c("80", "95"))
})

# ---------------------------------------------------------------------------
# Numeric sanity
# ---------------------------------------------------------------------------

test_that("predict.bets() mean values are finite", {
  fc <- predict(fit_ann, h = 5)
  expect_true(all(is.finite(fc$mean)))
})

test_that("predict.bets() x matches training series", {
  fc <- predict(fit_ann, h = 5)
  expect_identical(fc$x, fit_ann$y)
})

# ---------------------------------------------------------------------------
# Printing
# ---------------------------------------------------------------------------

test_that("print.bets_forecast() prints a table of forecasts and intervals", {
  out <- capture.output(print(predict(fit_seas, h = 4)))
  expect_length(out, 5)
  expect_match(out[1], "Point Forecast +Lo 80 +Hi 80 +Lo 95 +Hi 95")
  expect_match(out[2], "^7 Q1 ")
})

test_that("time_labels() gives months, quarters and years", {
  expect_identical(BETS:::time_labels(ts(1:3, start = c(2020, 11), frequency = 12)),
                   c("Nov 2020", "Dec 2020", "Jan 2021"))
  expect_identical(BETS:::time_labels(ts(1:2, start = c(2020, 4), frequency = 4)),
                   c("2020 Q4", "2021 Q1"))
  expect_identical(BETS:::time_labels(ts(1:2, start = 1999)), c("1999", "2000"))
})
