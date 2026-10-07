# ---------------------------------------------------------------------------
# Numerical robustness of the marginal likelihood
# ---------------------------------------------------------------------------

test_that("co2 ETS(A,A,A): no spurious evidence from exploding recursions", {
  # The recursion explodes for many theta; the old clamp of the quadratic form
  # at 0 gave log evidence ~ +600 instead of ~ -141
  mc <- BETS:::coerce_model_components("AAA", 12)
  ctrl <- BETS:::resolve_bets_control(list(N_final = 50L))
  le <- vapply(c("quadrature", "ais"), function(im) {
    set.seed(1)
    BETS:::fit_bets_models(co2, mc, ctrl, integration = im)$results[[1]]$log_evidence
  }, numeric(1))
  expect_true(all(le < -100))
  expect_lt(abs(le[["quadrature"]] - le[["ais"]]), 0.05)
})

test_that("lost-precision evaluations are invalid, not spikes (uncentred co2)", {
  # Without centring, no point may exceed the values of the centred series
  mc <- c("A", "A", "A", "FALSE")
  tn <- c("alpha", "beta", "gamma")
  make_lg <- function(y) {
    ctrl <- BETS:::resolve_bets_control(list())
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
  # From the heuristic point, the mode search stops 15.8 log units too low
  y <- nottem - mean(nottem[1:12])
  mc <- c("A", "A", "N", "FALSE")
  tn <- c("alpha", "beta")
  ctrl <- BETS:::resolve_bets_control(list())
  ctrl$psi0 <- 0.5 * (mean(diff(nottem, lag = 12)^2) + mean(diff(nottem)^2))
  lg <- BETS:::make_log_g_rb(y, mc, tn, ctrl, BETS:::init_rb_prior(y, mc, tn, ctrl))
  set.seed(1)
  best_of_1024 <- BETS:::prior_scan(lg, tn, n_scan = 1024)$log_g[1]
  set.seed(2)
  lap <- BETS:::laplace_mode(lg, BETS:::prior_scan(lg, tn, n_scan = 64)$Z[1, ])
  expect_gt(lg(matrix(lap$zhat, 1))$log_g, best_of_1024 - 1)
})

# ---------------------------------------------------------------------------
# Priors of sigma^2 and of the initial states
# ---------------------------------------------------------------------------

test_that("default_psi0() is the smaller residual variance of naive and seasonal naive", {
  set.seed(1)
  y <- ts(cumsum(rnorm(48)) + rep(10 * sin(2 * pi * (1:12) / 12), 4), frequency = 12)
  expect_equal(BETS:::default_psi0(y), min(var(diff(y)), var(diff(y, lag = 12))))
  # A drift does not count: variance, not mean square, of the differences
  y1 <- ts(cumsum(rnorm(30)) + 5 * (1:30))
  expect_equal(BETS:::default_psi0(y1), var(diff(y1)))
  # Positive for deterministic series
  expect_gt(BETS:::default_psi0(ts(1:20)), 0)
})

test_that("the initial-state prior depends on nu0 and psi0 only through E[sigma^2]", {
  set.seed(1)
  y <- ts(cumsum(rnorm(30)))
  mc <- c("A", "A", "N", "FALSE")
  tn <- c("alpha", "beta")
  V0 <- function(nu0, psi0) {
    ctrl <- BETS:::resolve_bets_control(list(nu0 = nu0, psi0 = psi0))
    BETS:::init_rb_prior(y, mc, tn, ctrl)$V0
  }
  expect_equal(V0(nu0 = 3, psi0 = 2), V0(nu0 = 10, psi0 = 16))   # E[sigma^2] = 2
  expect_equal(V0(nu0 = 3, psi0 = 4), V0(nu0 = 3, psi0 = 2) / 2)
})

test_that("bets() fits deterministic series (straight line, trend + exact seasonality)", {
  for (y in list(ts(1:20), ts(1:24 + rep(c(1, -1, 2, -2), 6), frequency = 4))) {
    set.seed(1)
    fit <- suppressWarnings(bets(y))
    expect_equal(sum(fit$fit$model_weights), 1, tolerance = 1e-8)
    expect_true(all(is.finite(predict(fit, h = 3)$mean)))
  }
})

test_that("centring the series leaves the evidence unchanged and restores the level", {
  set.seed(2)
  y <- ts(1000 + cumsum(rnorm(40)), frequency = 1)
  mc <- BETS:::coerce_model_components("AAN", 1)
  ctrl <- BETS:::resolve_bets_control(list())
  ctrl$psi0 <- mean(diff(y)^2)
  le_raw <- BETS:::quadrature_rb(y, mc[[1]], ctrl)$log_evidence
  set.seed(1)
  fit <- BETS:::fit_bets_models(y, mc, ctrl, integration = "quadrature")$results[[1]]
  # Equal up to the optimiser's tolerance on the mode (~1e-5).
  expect_lt(abs(fit$log_evidence - le_raw), 1e-3)
  # Final level states are on the scale of the data, not of the centred series.
  expect_lt(abs(mean(fit$states[, "l"]) - mean(tail(y, 5))), 5)
})
