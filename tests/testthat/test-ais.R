# ---------------------------------------------------------------------------
# adaptive_importance_sampling(): normalising constants of known targets
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

test_that("AIS recovers the normalising constant of a Gaussian (d = 1..4)", {
  # RQMC error with 1024 draws: at most 0.005 over 20 seeds for d = 4.
  set.seed(11)
  log_c <- -123.456
  for (d in 1:4) {
    A <- matrix(rnorm(d * d), d, d)
    S <- crossprod(A) / d + diag(0.2, d)
    mu <- seq(0.8, -0.6, length.out = d)
    ais <- BETS:::adaptive_importance_sampling(gaussian_log_g(mu, S, log_c),
                                               z_start = rep(0, d), n_draw = 1024,
                                               min_ess = 256)
    expect_equal(ais$n_iter, 1L)
    expect_lt(abs(ais$log_evidence - log_c), 0.02)
    expect_lt(max(abs(colSums(ais$w * ais$Z) - mu)), 0.05)
  }
})

test_that("AIS adapts and pools draws until min_ess is reached", {
  # Skewed target (Gumbel x logistic): one step of 256 draws has ESS of about
  # 200-230, below min_ess = 0.9 * 256, so a second iteration is needed; the
  # pooled estimate must stay unbiased.
  log_c <- -50
  skew_log_g <- function(Z) {
    list(log_g = log_c - (Z[, 1] + exp(-Z[, 1])) + stats::dlogis(Z[, 2], log = TRUE))
  }
  set.seed(1)
  ais <- BETS:::adaptive_importance_sampling(skew_log_g, z_start = c(0, 0),
                                             n_draw = 256, min_ess = 0.9 * 256)
  expect_gt(ais$n_iter, 1L)
  expect_gte(ais$ess, 0.9 * 256)
  expect_equal(nrow(ais$Z), ais$n_iter * 256)
  expect_lt(abs(ais$log_evidence - log_c), 0.02)
})

# ---------------------------------------------------------------------------
# adaptive_is_rb(): interface
# ---------------------------------------------------------------------------

test_that("adaptive_is_rb() returns posterior draws, ESS and proposal", {
  set.seed(3)
  y <- ts(cumsum(rnorm(30)) + 10, frequency = 4)
  for (mc in list(c("A", "N", "N", "FALSE"), c("A", "A", "A", "TRUE"))) {
    ctrl <- BETS:::resolve_bets_control(list(), 4)
    ctrl$psi0 <- mean(diff(y)^2)
    res <- BETS:::adaptive_is_rb(y, mc, ctrl)
    d <- if (mc[[2]] == "N") 1 else 4
    expect_named(res, c("thetas", "etas", "states", "sigma2s", "ess", "n_iter",
                        "prop_params", "log_evidence", "log_lik_pointwise", "timing"))
    expect_true(is.finite(res$log_evidence))
    expect_gte(res$ess, ctrl$N_draw[d] / 4)
    expect_gte(res$n_iter, 1L)
    expect_named(res$prop_params, c("mus", "Sigma", "df"))
    expect_equal(dim(res$thetas), c(ctrl$N_final[d], d))
  }
})

test_that("adaptive_is_rb() gives zero weight when min_ess is not reached", {
  set.seed(1)
  y <- ts(rnorm(20))
  ctrl <- BETS:::resolve_bets_control(list(min_ess = 1e6, N_iter_max = 2))
  ctrl$psi0 <- mean(diff(y)^2)
  expect_warning(res <- BETS:::adaptive_is_rb(y, c("A", "N", "N", "FALSE"), ctrl),
                 "zero weight")
  expect_equal(res$log_evidence, -Inf)
  expect_null(res$thetas)
  expect_equal(res$n_iter, 2L)
})
