# Helper: build a single-row theta_unc matrix with named columns
make_theta_unc <- function(vals, nms) {
  m <- matrix(vals, nrow = 1)
  colnames(m) <- nms
  m
}

# ---------------------------------------------------------------------------
# Alpha
# ---------------------------------------------------------------------------

test_that("alpha is in (0, 1) for negative unconstrained value", {
  tu <- make_theta_unc(-3, "alpha")
  res <- transform_unconstrained_to_theta(tu, "alpha", 0.8, 0.98)
  expect_true(res$theta[, "alpha"] > 0 && res$theta[, "alpha"] < 1)
})

test_that("alpha is in (0, 1) for zero unconstrained value", {
  tu <- make_theta_unc(0, "alpha")
  res <- transform_unconstrained_to_theta(tu, "alpha", 0.8, 0.98)
  expect_equal(as.numeric(res$theta[, "alpha"]), 0.5)
})

test_that("alpha is in (0, 1) for positive unconstrained value", {
  tu <- make_theta_unc(5, "alpha")
  res <- transform_unconstrained_to_theta(tu, "alpha", 0.8, 0.98)
  expect_true(res$theta[, "alpha"] > 0 && res$theta[, "alpha"] < 1)
})

test_that("alpha log-Jacobian is finite", {
  tu <- make_theta_unc(0, "alpha")
  res <- transform_unconstrained_to_theta(tu, "alpha", 0.8, 0.98)
  expect_true(is.finite(res$log_jac))
})

# ---------------------------------------------------------------------------
# Beta (depends on alpha: beta < alpha)
# ---------------------------------------------------------------------------

test_that("beta is strictly less than alpha", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "beta"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "beta"), 0.8, 0.98)
  expect_true(res$theta[, "beta"] < res$theta[, "alpha"])
})

test_that("beta is positive", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "beta"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "beta"), 0.8, 0.98)
  expect_true(res$theta[, "beta"] > 0)
})

test_that("beta log-Jacobian is finite", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "beta"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "beta"), 0.8, 0.98)
  expect_true(is.finite(res$log_jac))
})

# ---------------------------------------------------------------------------
# Gamma (depends on alpha: gamma < 1 - alpha)
# ---------------------------------------------------------------------------

test_that("gamma is strictly less than 1 - alpha", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "gamma"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "gamma"), 0.8, 0.98)
  alpha <- res$theta[, "alpha"]
  gamma <- res$theta[, "gamma"]
  expect_true(gamma < 1 - alpha)
})

test_that("gamma is positive", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "gamma"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "gamma"), 0.8, 0.98)
  expect_true(res$theta[, "gamma"] > 0)
})

test_that("gamma log-Jacobian is finite", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "gamma"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "gamma"), 0.8, 0.98)
  expect_true(is.finite(res$log_jac))
})

# ---------------------------------------------------------------------------
# Phi (bounded in [phi_min, phi_max])
# ---------------------------------------------------------------------------

test_that("phi is within (phi_min, phi_max) for unconstrained = 0", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "phi"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "phi"), 0.8, 0.98)
  phi <- res$theta[, "phi"]
  expect_true(phi > 0.8 && phi < 0.98)
})

test_that("phi converges to phi_max as unconstrained -> +Inf", {
  tu <- make_theta_unc(c(0, 100), c("alpha", "phi"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "phi"), 0.8, 0.98)
  expect_equal(as.numeric(res$theta[, "phi"]), 0.98, tolerance = 1e-6)
})

test_that("phi converges to phi_min as unconstrained -> -Inf", {
  tu <- make_theta_unc(c(0, -100), c("alpha", "phi"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "phi"), 0.8, 0.98)
  expect_equal(as.numeric(res$theta[, "phi"]), 0.8, tolerance = 1e-6)
})

test_that("phi log-Jacobian is finite", {
  tu <- make_theta_unc(c(0, 0), c("alpha", "phi"))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "phi"), 0.8, 0.98)
  expect_true(is.finite(res$log_jac))
})

# ---------------------------------------------------------------------------
# Full 4-parameter transform
# ---------------------------------------------------------------------------

test_that("full 4-param transform returns finite log_jac", {
  tu <- make_theta_unc(c(0, 0, 0, 0), c("alpha", "beta", "gamma", "phi"))
  res <- transform_unconstrained_to_theta(
    tu, c("alpha", "beta", "gamma", "phi"), 0.8, 0.98
  )
  expect_true(is.finite(res$log_jac))
})

test_that("full 4-param transform satisfies all constraints simultaneously", {
  tu <- make_theta_unc(c(0.5, -0.3, 0.2, 1.0), c("alpha", "beta", "gamma", "phi"))
  res <- transform_unconstrained_to_theta(
    tu, c("alpha", "beta", "gamma", "phi"), 0.8, 0.98
  )
  theta <- res$theta
  alpha <- theta[, "alpha"]
  beta  <- theta[, "beta"]
  gamma <- theta[, "gamma"]
  phi   <- theta[, "phi"]
  expect_true(alpha > 0 && alpha < 1)
  expect_true(beta  > 0 && beta  < alpha)
  expect_true(gamma > 0 && gamma < 1 - alpha)
  expect_true(phi   > 0.8 && phi < 0.98)
})

test_that("transform works row-wise for multiple particles", {
  set.seed(7)
  N <- 10
  tu <- matrix(rnorm(N * 2), nrow = N, ncol = 2,
               dimnames = list(NULL, c("alpha", "beta")))
  res <- transform_unconstrained_to_theta(tu, c("alpha", "beta"), 0.8, 0.98)
  expect_equal(nrow(res$theta), N)
  expect_length(res$log_jac, N)
  expect_true(all(res$theta[, "beta"] < res$theta[, "alpha"]))
  expect_true(all(is.finite(res$log_jac)))
})

test_that("transform_theta_to_unconstrained() inverts the transform", {
  set.seed(8)
  pars <- c("alpha", "beta", "phi", "gamma")
  z <- matrix(rnorm(40, sd = 3), ncol = 4, dimnames = list(NULL, pars))
  theta <- transform_unconstrained_to_theta(z, pars, 0.8, 0.98)$theta
  expect_equal(transform_theta_to_unconstrained(theta, 0.8, 0.98), z, tolerance = 1e-8)
  expect_equal(transform_theta_to_unconstrained(theta[, "alpha", drop = FALSE], 0.8, 0.98),
               z[, "alpha", drop = FALSE], tolerance = 1e-8)
})
