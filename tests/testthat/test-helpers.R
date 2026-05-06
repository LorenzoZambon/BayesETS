# ---------------------------------------------------------------------------
# ets_label()
# ---------------------------------------------------------------------------

test_that("ets_label() returns 'ANN' for ANN model", {
  expect_equal(BETS:::ets_label(c("A", "N", "N", "FALSE")), "ANN")
})

test_that("ets_label() returns 'AAN' for trend model", {
  expect_equal(BETS:::ets_label(c("A", "A", "N", "FALSE")), "AAN")
})

test_that("ets_label() returns 'AAdN' for damped trend model", {
  expect_equal(BETS:::ets_label(c("A", "A", "N", "TRUE")), "AAdN")
})

test_that("ets_label() returns 'AAA' for full additive model", {
  expect_equal(BETS:::ets_label(c("A", "A", "A", "FALSE")), "AAA")
})

test_that("ets_label() returns 'AAdA' for damped seasonal model", {
  expect_equal(BETS:::ets_label(c("A", "A", "A", "TRUE")), "AAdA")
})

test_that("ets_label() returns 'ANA' for seasonal-only model", {
  expect_equal(BETS:::ets_label(c("A", "N", "A", "FALSE")), "ANA")
})

# ---------------------------------------------------------------------------
# compute_stacking_weights() — validation
# ---------------------------------------------------------------------------

test_that("compute_stacking_weights() errors on empty list", {
  expect_error(BETS:::compute_stacking_weights(list()), "log_lik_list is empty")
})

test_that("compute_stacking_weights() errors when element is not a matrix", {
  expect_error(
    BETS:::compute_stacking_weights(list(1:10)),
    "must be a matrix"
  )
})

test_that("compute_stacking_weights() errors on mismatched column counts", {
  expect_error(
    BETS:::compute_stacking_weights(list(
      matrix(0, nrow = 5, ncol = 10),
      matrix(0, nrow = 5, ncol = 8)
    )),
    "same number of columns"
  )
})

# ---------------------------------------------------------------------------
# compute_stacking_weights() — correctness
# ---------------------------------------------------------------------------

test_that("compute_stacking_weights() returns weights summing to 1", {
  set.seed(1)
  ll1 <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  ll2 <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
  w <- BETS:::compute_stacking_weights(list(ll1, ll2))
  expect_length(w, 2)
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

test_that("compute_stacking_weights() returns non-negative weights", {
  set.seed(2)
  ll1 <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  ll2 <- matrix(rnorm(50 * 15), nrow = 50, ncol = 15)
  w <- BETS:::compute_stacking_weights(list(ll1, ll2))
  expect_true(all(w >= 0))
})

test_that("compute_stacking_weights() works with a single model (weight = 1)", {
  ll1 <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  w <- BETS:::compute_stacking_weights(list(ll1))
  expect_equal(w, 1, tolerance = 1e-6)
})

test_that("compute_stacking_weights() handles K=3 models", {
  set.seed(3)
  mats <- lapply(1:3, function(i) matrix(rnorm(80 * 20), nrow = 80, ncol = 20))
  w <- BETS:::compute_stacking_weights(mats)
  expect_length(w, 3)
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# simulate_future_trajectories()
# ---------------------------------------------------------------------------

# Shared fit object used across trajectory tests
local({
  set.seed(42)
  fit_ann <<- bets(ts(rnorm(30)), model = "ANN")
})

test_that("simulate_future_trajectories() returns a matrix", {
  traj <- BETS:::simulate_future_trajectories(fit_ann$fit, h = 5, n_traj = 50)
  expect_true(is.matrix(traj))
})

test_that("simulate_future_trajectories() returns n_traj rows", {
  traj <- BETS:::simulate_future_trajectories(fit_ann$fit, h = 5, n_traj = 50)
  expect_equal(nrow(traj), 50)
})

test_that("simulate_future_trajectories() returns h columns", {
  traj <- BETS:::simulate_future_trajectories(fit_ann$fit, h = 8, n_traj = 40)
  expect_equal(ncol(traj), 8)
})

test_that("simulate_future_trajectories() returns finite values", {
  traj <- BETS:::simulate_future_trajectories(fit_ann$fit, h = 5, n_traj = 50)
  expect_true(all(is.finite(traj)))
})
