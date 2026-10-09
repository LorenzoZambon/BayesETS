# ---------------------------------------------------------------------------
# fitted() and residuals()
# ---------------------------------------------------------------------------

local({
  set.seed(1)
  y <- ts(10 + cumsum(rnorm(40)) + rep(c(2, -1, 0, -1), 10), start = c(2000, 1), frequency = 4)
  fit_fr <<- bets(y)
})

test_that("type = 'mean': ts with the time index of y, fitted + residuals = y", {
  f <- fitted(fit_fr)
  r <- residuals(fit_fr)
  expect_true(stats::is.ts(f))
  expect_true(stats::is.ts(r))
  expect_equal(stats::tsp(f), stats::tsp(fit_fr$y))
  expect_equal(as.numeric(f + r), as.numeric(fit_fr$y))
})

test_that("type = 'draws': n_draws rows, consistent with the mean", {
  set.seed(2)
  rd <- residuals(fit_fr, type = "draws", n_draws = 4000)
  set.seed(2)
  fd <- fitted(fit_fr, type = "draws", n_draws = 4000)
  L <- length(fit_fr$y)
  expect_equal(dim(rd), c(4000L, L))
  expect_equal(fd + rd, matrix(as.numeric(fit_fr$y), 4000, L, byrow = TRUE))
  expect_lt(max(abs(colMeans(rd) - as.numeric(residuals(fit_fr)))), 0.1 * sd(fit_fr$y))
})

test_that("the mean uses the model weights", {
  set.seed(3)
  y <- ts(cumsum(rnorm(50)))
  set.seed(1)
  fit_one <- bets(y, model = "ANN")
  set.seed(1)
  fit_two <- bets(y, model = list("ANN", "AAN"), control = list(prior_models = c(1, 0)))
  expect_equal(residuals(fit_two), residuals(fit_one))
  # One-step residuals of a random walk are close to its differences
  expect_gt(cor(as.numeric(residuals(fit_one))[-1], diff(as.numeric(y))), 0.9)
})

test_that("constant series: zero residuals", {
  fit <- suppressWarnings(bets(ts(rep(5, 12))))
  expect_equal(as.numeric(residuals(fit)), rep(0, 12))
  expect_equal(as.numeric(fitted(fit)), rep(5, 12))
})

test_that("fitted() and residuals() check their arguments", {
  expect_error(residuals(fit_fr, type = "median"), "should be one of")
  expect_error(fitted(fit_fr, type = "draws", n_draws = 0), "n_draws must be")
  expect_error(fitted(fit_fr, type = "draws", n_draws = 10.5), "n_draws must be")
})
