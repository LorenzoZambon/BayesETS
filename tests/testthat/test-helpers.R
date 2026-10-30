# ---------------------------------------------------------------------------
# ets_label()
# ---------------------------------------------------------------------------

test_that("ets_label() returns 'ANN' for ANN model", {
  expect_equal(ets_label(c("A", "N", "N", "FALSE")), "ANN")
})

test_that("ets_label() returns 'AAN' for trend model", {
  expect_equal(ets_label(c("A", "A", "N", "FALSE")), "AAN")
})

test_that("ets_label() returns 'AAdN' for damped trend model", {
  expect_equal(ets_label(c("A", "A", "N", "TRUE")), "AAdN")
})

test_that("ets_label() returns 'AAA' for full additive model", {
  expect_equal(ets_label(c("A", "A", "A", "FALSE")), "AAA")
})

test_that("ets_label() returns 'AAdA' for damped seasonal model", {
  expect_equal(ets_label(c("A", "A", "A", "TRUE")), "AAdA")
})

test_that("ets_label() returns 'ANA' for seasonal-only model", {
  expect_equal(ets_label(c("A", "N", "A", "FALSE")), "ANA")
})

# ---------------------------------------------------------------------------
# compute_stacking_weights(): validation
# ---------------------------------------------------------------------------

test_that("compute_stacking_weights() errors on empty list", {
  expect_error(compute_stacking_weights(list()), "log_lik_list is empty")
})

test_that("compute_stacking_weights() errors when element is not a matrix", {
  expect_error(
    compute_stacking_weights(list(1:10)),
    "must be a matrix"
  )
})

test_that("compute_stacking_weights() errors on mismatched column counts", {
  expect_error(
    compute_stacking_weights(list(
      matrix(0, nrow = 5, ncol = 10),
      matrix(0, nrow = 5, ncol = 8)
    )),
    "same number of columns"
  )
})

# ---------------------------------------------------------------------------
# compute_stacking_weights(): correctness
# ---------------------------------------------------------------------------

test_that("compute_stacking_weights() returns weights summing to 1", {
  set.seed(1)
  ll1 <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  ll2 <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  w <- compute_stacking_weights(list(ll1, ll2))
  expect_length(w, 2)
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("compute_stacking_weights() returns non-negative weights", {
  set.seed(2)
  ll1 <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  ll2 <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  w <- compute_stacking_weights(list(ll1, ll2))
  expect_true(all(w >= 0))
})

test_that("compute_stacking_weights() works with a single model (weight = 1)", {
  ll1 <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  w <- compute_stacking_weights(list(ll1))
  expect_equal(w, 1, tolerance = 1e-6)
})

test_that("compute_stacking_weights() handles K=3 models", {
  set.seed(3)
  mats <- lapply(1:3, function(i) matrix(rnorm(80 * 20), nrow = 80, ncol = 20))
  w <- compute_stacking_weights(mats)
  expect_length(w, 3)
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# simulate_paths()
# ---------------------------------------------------------------------------

# Shared fit object used across path tests
local({
  set.seed(42)
  fit_ann <<- bets(ts(rnorm(30)), model = "ANN")
})

test_that("simulate_paths() returns a finite n_traj x h matrix", {
  for (h in c(1L, 8L)) {
    paths <- simulate_paths(fit_ann$fit, forecast_mixture(fit_ann$fit, h), h, 50L)
    expect_true(is.matrix(paths))
    expect_equal(dim(paths), c(50L, h))
    expect_true(all(is.finite(paths)))
  }
  expect_equal(dim(simulate_paths(fit_ann$fit, forecast_mixture(fit_ann$fit, 3L), 3L, 1L)), c(1L, 3L))
})

# ---------------------------------------------------------------------------
# heavy_set()
# ---------------------------------------------------------------------------

test_that("heavy_set() keeps the smallest set of largest weights with the mass", {
  w <- c(0.5, 1e-6, 0.3, 2e-6, 0.2 - 3e-6)
  expect_identical(heavy_set(w, 1e-5), c(1L, 3L, 5L))
  expect_identical(heavy_set(w, 0), 1:5)
  expect_identical(heavy_set(w, 0.25), c(1L, 3L))
  expect_identical(heavy_set(c(0, 1, 0), 1e-4), 2L)
  set.seed(1)
  w <- stats::rexp(3000)^4
  k <- heavy_set(w, 1e-4)
  expect_gte(sum(w[k]), (1 - 1e-4) * sum(w))
  expect_lt(sum(w[k]) - min(w[k]), (1 - 1e-4) * sum(w))   # minimal
  expect_false(is.unsorted(k))
})
