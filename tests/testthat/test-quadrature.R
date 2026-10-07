# ---------------------------------------------------------------------------
# gauss_hermite_prob(): probabilists' convention
# ---------------------------------------------------------------------------

test_that("converted Gauss-Hermite weights sum to 1", {
  for (k in c(3, 5, 7, 9)) {
    gh <- BETS:::gauss_hermite_prob(k)
    expect_equal(sum(gh$w), 1, tolerance = 1e-12)
  }
})

test_that("converted Gauss-Hermite rule integrates N(0, 1) moments", {
  # Catches a missing sqrt(2) on the nodes: E[X^2] = 1, E[X^4] = 3.
  for (k in c(3, 5, 7, 9)) {
    gh <- BETS:::gauss_hermite_prob(k)
    expect_equal(sum(gh$w * gh$x^2), 1, tolerance = 1e-12)
    expect_equal(sum(gh$w * gh$x^4), 3, tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# adaptive_gh_quadrature(): normalising constant of a known Gaussian
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

test_that("quadrature recovers the normalising constant of a Gaussian (d = 1..4)", {
  set.seed(11)
  log_c <- -123.456
  for (d in 1:4) {
    A <- matrix(rnorm(d * d), d, d)
    S <- crossprod(A) / d + diag(0.2, d)
    mu <- seq(0.8, -0.6, length.out = d)
    q <- BETS:::adaptive_gh_quadrature(gaussian_log_g(mu, S, log_c),
                                       z_start = rep(0, d), n_quad = 5)
    expect_lt(abs(q$log_evidence - log_c), 1e-8)
    # Normalised weights reproduce the mean of the Gaussian.
    expect_lt(max(abs(colSums(q$w * q$Z) - mu)), 1e-8)
  }
})

test_that("quadrature keeps matrices as matrices for d = 1", {
  q <- BETS:::adaptive_gh_quadrature(gaussian_log_g(0.3, matrix(0.5), 0),
                                     z_start = 0, n_quad = 7)
  expect_true(is.matrix(q$Z))
  expect_equal(dim(q$Z), c(7L, 1L))
  expect_equal(dim(q$H), c(1L, 1L))
})

test_that("quadrature warns and floors a non positive-definite Hessian", {
  # Flat direction in z2: -log g has zero curvature there.
  log_g_flat <- function(Z) list(log_g = -0.5 * Z[, 1]^2)
  expect_warning(
    q <- BETS:::adaptive_gh_quadrature(log_g_flat, z_start = c(0, 0), n_quad = 3),
    "not safely positive definite"
  )
  expect_true(is.finite(q$log_evidence))
})

# ---------------------------------------------------------------------------
# fit_one_model() with quadrature: interface
# ---------------------------------------------------------------------------

test_that("fit_one_model() returns the same fields with quadrature and AIS", {
  set.seed(3)
  y <- ts(cumsum(rnorm(30)) + 10, frequency = 4)
  for (mc in list(c("A", "N", "N", "FALSE"), c("A", "A", "A", "TRUE"))) {
    ctrl <- BETS:::resolve_bets_control(list())
    ctrl$psi0 <- mean(diff(y)^2)
    res_ais  <- BETS:::fit_one_model(y, mc, ctrl, "ais")
    res_quad <- BETS:::fit_one_model(y, mc, ctrl, "quadrature")
    expect_named(res_quad, names(res_ais))
    expect_true(is.na(res_quad$ess))
    expect_true(is.na(res_quad$n_iter))
    expect_null(res_quad$prop_params)
    expect_true(is.finite(res_quad$log_evidence))
    expect_equal(dim(res_quad$thetas), dim(res_ais$thetas))
    expect_equal(colnames(res_quad$thetas), colnames(res_ais$thetas))
    expect_equal(colnames(res_quad$states), colnames(res_ais$states))
  }
})

test_that("bets() with control integration = 'quadrature' fits and predicts", {
  set.seed(1)
  fit <- bets(ts(rnorm(20)), model = "AAdN", control = list(integration = "quadrature"))
  expect_equal(fit$fit$results[[1]]$integration, "quadrature")
  expect_true(is.na(fit$fit$results[[1]]$ess))
  fc <- predict(fit, h = 5)
  expect_true(all(is.finite(fc$mean)))
  expect_true(all(fc$lower < fc$upper))
})

test_that("quadrature works with stacking (pointwise log-likelihoods)", {
  set.seed(1)
  y <- ts(rnorm(30))
  ctrl <- BETS:::resolve_bets_control(list(N_final = 50L))
  fit <- BETS:::fit_bets_models(y, BETS:::coerce_model_components(list("ANN", "AAN"), 1),
                                ctrl, combination = "stacking", integration = "quadrature")
  expect_equal(dim(fit$results[[1]]$log_lik_pointwise), c(50L, 30L))
  expect_equal(sum(fit$model_weights), 1, tolerance = 1e-8)
})
