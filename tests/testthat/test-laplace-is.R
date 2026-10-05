# ---------------------------------------------------------------------------
# laplace_importance_sampling(): normalising constant of a known Gaussian
# ---------------------------------------------------------------------------

# log g(z) = log_c + log N(z; mu, S), so that integral g(z) dz = exp(log_c).
gaussian_log_g <- function(mu, S, log_c) {
  R <- chol(S)
  d <- length(mu)
  function(Z) {
    dev <- sweep(Z, 2, mu, "-")
    mahal <- colSums(backsolve(R, t(dev), transpose = TRUE)^2)
    list(log_g = log_c - 0.5 * mahal - sum(log(diag(R))) - (d / 2) * log(2 * pi))
  }
}

test_that("Laplace IS recovers the normalising constant of a Gaussian (d = 1..4)", {
  # RQMC error with 1024 draws: at most 0.005 over 20 seeds for d = 4.
  set.seed(11)
  log_c <- -123.456
  for (d in 1:4) {
    A <- matrix(rnorm(d * d), d, d)
    S <- crossprod(A) / d + diag(0.2, d)
    mu <- seq(0.8, -0.6, length.out = d)
    lis <- BETS:::laplace_importance_sampling(gaussian_log_g(mu, S, log_c),
                                              z_start = rep(0, d), n_is = 1024)
    expect_lt(abs(lis$log_evidence - log_c), 0.02)
    expect_lt(max(abs(colSums(lis$w * lis$Z) - mu)), 0.05)
    expect_gt(lis$ess, 0.5 * 1024)
  }
})

# ---------------------------------------------------------------------------
# laplace_is_rb(): interface
# ---------------------------------------------------------------------------

test_that("laplace_is_rb() returns the same fields as adaptive_is_rb()", {
  set.seed(3)
  y <- ts(cumsum(rnorm(30)) + 10, frequency = 4)
  for (mc in list(c("A", "N", "N", "FALSE"), c("A", "A", "A", "TRUE"))) {
    ctrl <- BETS:::resolve_bets_control(list(), 4)
    ctrl$psi0 <- mean(diff(y)^2)
    res_ais <- BETS:::adaptive_is_rb(y, mc, ctrl)
    res_lis <- BETS:::laplace_is_rb(y, mc, ctrl)
    expect_named(res_lis, names(res_ais))
    expect_true(is.finite(res_lis$log_evidence))
    expect_true(res_lis$ess > 0)
    expect_true(is.na(res_lis$n_iter))
    expect_named(res_lis$prop_params, names(res_ais$prop_params))
    expect_equal(dim(res_lis$thetas), dim(res_ais$thetas))
    expect_equal(colnames(res_lis$thetas), colnames(res_ais$thetas))
    expect_equal(colnames(res_lis$states), colnames(res_ais$states))
  }
})

test_that("bets() dispatches to laplace_is and predicts", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "ANN",
              control = list(integration_method = "laplace_is"))
  expect_true(is.finite(fit$fit$results[[1]]$log_evidence))
  fc <- predict(fit, h = 5)
  expect_true(all(is.finite(fc$mean)))
  expect_true(all(fc$lower < fc$upper))
})
