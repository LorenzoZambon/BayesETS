# ---------------------------------------------------------------------------
# Numerical robustness of the marginal likelihood
# ---------------------------------------------------------------------------

test_that("co2 ETS(A,A,A): no spurious evidence from exploding recursions", {
  # Long, very smooth series: many theta in the prior support make the
  # recursion explode.  The residual quadratic form used to be clamped at 0
  # there, which gave log evidence ~ +600 instead of ~ -196.
  mc <- BETS:::coerce_model_components("AAA", 12)
  ctrl <- BETS:::resolve_bets_control(list(N_final = 50L), 12)
  le <- vapply(c("quadrature", "ais"), function(im) {
    set.seed(1)
    BETS:::fit_bets_models(co2, mc, ctrl, integration = im)$results[[1]]$log_evidence
  }, numeric(1))
  expect_true(all(le < -150))
  expect_lt(abs(le[["quadrature"]] - le[["ais"]]), 0.05)
})

test_that("lost-precision evaluations are invalid, not spikes (uncentred co2)", {
  # At unstable theta the residual quadratic form loses all its digits; the
  # old clamp at 0 turned these points into log g ~ +600.  Without centring,
  # no point may exceed the values of the (numerically safe) centred series.
  mc <- c("A", "A", "A", "FALSE")
  tn <- c("alpha", "beta", "gamma")
  make_lg <- function(y) {
    ctrl <- BETS:::resolve_bets_control(list(), 12)
    ctrl$psi0 <- 0.5 * (mean(diff(y, lag = 12)^2) + mean(diff(y)^2))
    BETS:::make_log_g_rb(y, mc, tn, ctrl, BETS:::init_rb_prior(y, mc, tn, ctrl))
  }
  set.seed(1)
  Z <- matrix(runif(3 * 500, -3, 3), ncol = 3)
  lg_raw <- make_lg(co2)(Z)$log_g
  lg_centred <- make_lg(co2 - mean(co2[1:12]))(Z)$log_g
  expect_lt(max(lg_raw), max(lg_centred) + 1)
})

test_that("prior_scan() returns the candidates sorted by log g", {
  tn <- c("alpha", "beta", "gamma")
  log_g <- function(Z) list(log_g = -rowSums((Z - 0.3)^2))
  set.seed(1)
  s <- BETS:::prior_scan(log_g, tn, n_scan = 16)
  expect_equal(dim(s$Z), c(18L, 3L))   # heuristic + small smoothing + 16 Sobol points
  expect_false(is.unsorted(rev(s$log_g)))
  expect_true(any(apply(s$Z, 1, function(z) all(z == BETS:::heuristic_z_start(tn)))))
  expect_true(any(apply(s$Z, 1, function(z) all(z == BETS:::small_smoothing_z(tn)))))
})

test_that("the mode search is not trapped in the unstable region (nottem, AAN)", {
  # From the heuristic start, L-BFGS-B stops 15.8 log units below the best
  # point of a 1024-point prior scan; from the default 64-point scan it does not.
  y <- nottem - mean(nottem[1:12])
  mc <- c("A", "A", "N", "FALSE")
  tn <- c("alpha", "beta")
  ctrl <- BETS:::resolve_bets_control(list(), 12)
  ctrl$psi0 <- 0.5 * (mean(diff(nottem, lag = 12)^2) + mean(diff(nottem)^2))
  lg <- BETS:::make_log_g_rb(y, mc, tn, ctrl, BETS:::init_rb_prior(y, mc, tn, ctrl))
  set.seed(1)
  best_of_1024 <- BETS:::prior_scan(lg, tn, n_scan = 1024)$log_g[1]
  set.seed(2)
  lap <- BETS:::laplace_mode(lg, BETS:::prior_scan(lg, tn, n_scan = 64)$Z[1, ])
  expect_gt(lg(matrix(lap$zhat, 1))$log_g, best_of_1024 - 1)
})

test_that("centring the series leaves the evidence unchanged and restores the level", {
  set.seed(2)
  y <- ts(1000 + cumsum(rnorm(40)), frequency = 1)
  mc <- BETS:::coerce_model_components("AAN", 1)
  ctrl <- BETS:::resolve_bets_control(list(), 1)
  ctrl$psi0 <- mean(diff(y)^2)
  le_raw <- BETS:::quadrature_rb(y, mc[[1]], ctrl)$log_evidence
  set.seed(1)
  fit <- BETS:::fit_bets_models(y, mc, ctrl, integration = "quadrature")$results[[1]]
  # Equal up to the optimiser's tolerance on the mode (~1e-5).
  expect_lt(abs(fit$log_evidence - le_raw), 1e-3)
  # Final level states are on the scale of the data, not of the centred series.
  expect_lt(abs(mean(fit$states[, "l"]) - mean(tail(y, 5))), 5)
})
