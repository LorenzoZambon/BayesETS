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

test_that("predict.bets() default h: the frequency if in (1, 24], 10 otherwise", {
  expect_length(predict(fit_ann, n_traj = 50)$mean, 10)
  expect_length(predict(fit_seas, n_traj = 50)$mean, 4)
  set.seed(42)
  fit_48 <- bets(ts(rnorm(60), frequency = 48), model = "ANN")
  expect_length(predict(fit_48, n_traj = 50)$mean, 10)
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

test_that("predict.bets(return_paths = TRUE) also returns sample paths", {
  set.seed(1)
  fc <- predict(fit_seas, h = 6, level = c(50, 90), n_traj = 400, return_paths = TRUE)
  expect_true(stats::is.ts(fc$paths))
  expect_equal(dim(fc$paths), c(6L, 400L))
  expect_equal(stats::tsp(fc$paths), stats::tsp(fc$mean))
  # The forecasts do not depend on the paths
  fc0 <- predict(fit_seas, h = 6, level = c(50, 90))
  expect_identical(fc0[c("mean", "lower", "upper")], fc[c("mean", "lower", "upper")])
  # Default, h = 1, invalid flag
  expect_null(predict(fit_seas, h = 2)$paths)
  expect_equal(dim(predict(fit_ann, h = 1, n_traj = 10, return_paths = TRUE)$paths), c(1L, 10L))
  expect_error(predict(fit_ann, return_paths = NA), "return_paths must be TRUE or FALSE")
})

test_that("predict.bets(return_mixture = TRUE) returns the t mixture behind the forecasts", {
  fc <- predict(fit_seas, h = 6, level = c(80, 95), return_mixture = TRUE)
  mix <- fc$mixture
  expect_named(mix, c("weights", "df", "location", "scale"))
  K <- length(mix$weights)
  expect_equal(sum(mix$weights), 1)
  expect_equal(dim(mix$location), c(K, 6L))
  expect_equal(dim(mix$scale), c(K, 6L))
  expect_equal(mix$df, fit_seas$control$nu0 + length(fit_seas$y))
  expect_equal(as.numeric(fc$mean), colSums(mix$weights * mix$location))
  # The bounds solve the mixture CDF
  cdf <- function(x, j) sum(mix$weights * stats::pt((x - mix$location[, j]) / mix$scale[, j], mix$df))
  for (j in 1:6) {
    expect_equal(c(cdf(fc$lower[j, "95"], j), cdf(fc$lower[j, "80"], j), cdf(fc$upper[j, "80"], j),
                   cdf(fc$upper[j, "95"], j)), c(0.025, 0.1, 0.9, 0.975), tolerance = 1e-7)
  }
  # Both flags; default output without them; invalid flag
  fc2 <- predict(fit_seas, h = 6, return_paths = TRUE, return_mixture = TRUE, n_traj = 20)
  expect_true(all(c("paths", "mixture") %in% names(fc2)))
  expect_null(predict(fit_seas, h = 6)$mixture)
  expect_error(predict(fit_ann, return_mixture = 1), "return_mixture must be TRUE or FALSE")
})

test_that("predict.bets() is deterministic without paths and does not use random numbers", {
  set.seed(3)
  seed <- get(".Random.seed", envir = globalenv())
  fc1 <- predict(fit_seas, h = 9)
  fc2 <- predict(fit_seas, h = 9)
  expect_identical(fc1[c("mean", "lower", "upper")], fc2[c("mean", "lower", "upper")])
  expect_identical(get(".Random.seed", envir = globalenv()), seed)
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
  expect_identical(time_labels(ts(1:3, start = c(2020, 11), frequency = 12)),
                   c("Nov 2020", "Dec 2020", "Jan 2021"))
  expect_identical(time_labels(ts(1:2, start = c(2020, 4), frequency = 4)),
                   c("2020 Q4", "2021 Q1"))
  expect_identical(time_labels(ts(1:2, start = 1999)), c("1999", "2000"))
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("plot.bets_forecast() plots the series and the forecasts, and returns them invisibly", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  fc <- predict(fit_seas, h = 4)
  expect_false(withVisible(plot(fc))$visible)
  expect_identical(plot(fc), fc)
  # x range: from the first of the last `include` observations to the end of the forecasts
  plot(fc, include = 8)
  xlim <- c(time(fc$x)[length(fc$x) - 7], time(fc$mean)[4])
  expect_equal(graphics::par("usr")[1:2], xlim + c(-1, 1) * 0.04 * diff(xlim))
  # Graphical parameters in ... override the defaults
  plot(fc, ylim = c(-10, 10))
  expect_equal(graphics::par("usr")[3:4], c(-10.8, 10.8))
  expect_no_error(plot(predict(fit_ann, h = 1, level = 90)))
  expect_error(plot(fc, include = 0), "include must be")
})

test_that("interval_shade(): lighter for higher levels, between white and the colour", {
  shades <- interval_shade(c(10, 50, 80, 95, 99), "steelblue")
  brightness <- colSums(grDevices::col2rgb(shades))
  expect_true(all(diff(brightness) > 0))
  expect_true(all(brightness > sum(grDevices::col2rgb("steelblue")) & brightness < 3 * 255))
  expect_identical(interval_shade(95, "steelblue"), shades[4])
})

# ---------------------------------------------------------------------------
# Closed-form forecast distribution
# ---------------------------------------------------------------------------

test_that("final states: posterior mean and factor agree with the ETS recursion", {
  # Final states from the initial states eta, by ets_residuals(): at eta = mu_n they are
  # final_mean; the shift by row r of Rn moves them by row r of final_R
  set.seed(1)
  y <- ts(10 + rep(c(2, -1, 0, -1), 8) + cumsum(rnorm(32, sd = 0.3)), frequency = 4)
  for (code in c("ANN", "AAdN", "AAdA")) {
    fit <- bets(y, model = code)
    r <- fit$fit$results[[1]]
    p <- r$particles
    mc <- r$model_components
    flags <- model_flags(mc)
    m <- fit$period
    yc <- y - centring_shift(y)
    ctrl <- utils::modifyList(fit$control, list(psi0 = fit$psi0, s2 = fit$s2))
    tn <- colnames(p$theta)
    ml <- make_log_g_rb(yc, mc, tn, ctrl, init_rb_prior(yc, mc, tn, ctrl))(
      transform_theta_to_unconstrained(p$theta, ctrl$phi_min, ctrl$phi_max))$ml_res
    cpp_names <- c("l", if (flags$trend) "b", if (flags$seas) paste0("s", rev(seq_len(m))))
    perm <- match(state_names_of(mc, m), cpp_names)
    n <- length(perm)
    final <- function(eta) {
      t(ets_residuals(as.numeric(yc), flags$trend, flags$seas, flags$damped, m,
                      t(eta[perm, , drop = FALSE]), p$theta)$states)
    }
    x0 <- final(ml$mu_n)
    x0[1, ] <- x0[1, ] + centring_shift(y)
    expect_equal(unname(p$final_mean), unname(x0), tolerance = 1e-8, info = code)
    for (k in seq_len(n)) {
      d <- final(ml$mu_n + matrix(ml$Rn[k, , ], nrow = n)) - final(ml$mu_n)
      expect_equal(unname(matrix(p$final_R[k, , ], nrow = n)), unname(d), tolerance = 1e-8, info = code)
    }
    # h = 1: location of each component = one-step forecast from its final states
    one_step <- x0[1, ] + (if (flags$trend) p$theta[, "phi"] * x0[2, ] else 0) +
      (if (flags$seas) x0["s1" == state_names_of(mc, m), ] else 0)
    expect_equal(ets_t_mixture(r, 1L, seq_along(p$w))$location[, 1], unname(one_step), info = code)
  }
})

test_that("forecast mixture: weights of models and particles, truncation", {
  set.seed(1)
  y <- ts(10 + rep(c(2, -1, 0, -1), 8) + rnorm(32, sd = 0.5), frequency = 4)
  fit <- bets(y, model = c("ANN", "ANA", "AAdA"))
  mix <- forecast_mixture(fit$fit, 5L)
  w <- lapply(fit$fit$results, function(r) r$particles$w)
  W <- fit$fit$model_weights[mix$model] * mapply(function(i, k) w[[i]][k], mix$model, mix$particle)
  expect_equal(mix$weights, W / sum(W))
  expect_equal(dim(mix$location), c(length(mix$weights), 5L))
  # Models with final states: weight >= 1e-4; their components left out hold <= 0.01%
  with_final <- which(vapply(fit$fit$results, function(r) !is.null(r$particles$final_mean), logical(1)))
  expect_true(all(fit$fit$model_weights[with_final] >= 1e-4))
  expect_true(all(fit$fit$model_weights[-with_final] < 1e-4))
  expect_gte(sum(W), (1 - 1e-4) * sum(fit$fit$model_weights[with_final]))
})

test_that("negligible particles and models: dropped from the forecasts, not from the summaries", {
  # Deterministic quadrature fit (no random prior scan): the particles kept hold all
  # but at most 1e-5 of the weight of the nodes
  set.seed(1)
  y <- ts(cumsum(rnorm(40)))
  ctrl <- resolve_bets_control(list(n_scan = 0L))
  ctrl <- utils::modifyList(ctrl, prior_scales(y, ctrl))
  mc <- coerce_model_components("AAN", 1)[[1]]
  tn <- theta_names_of(mc)
  res <- fit_one_model(y, mc, ctrl, "quadrature")
  prior <- init_rb_prior(y, mc, tn, ctrl)
  log_g_value <- make_log_g_rb(y, mc, tn, ctrl, prior, posterior = FALSE)
  int <- integrate_quadrature(make_log_g_rb(y, mc, tn, ctrl, prior), prior_scan(log_g_value, tn, 0)$Z[1, ],
                              tn, ctrl, log_g_mode = log_g_value)
  keep <- heavy_set(int$w, 1e-5)
  expect_lt(length(keep), length(int$w))
  expect_lte(sum(int$w[-keep]), 1e-5)
  expect_equal(res$particles$theta, int$theta[keep, , drop = FALSE])
  # A model with zero weight has no final states, but keeps its posterior
  fit <- bets(ts(rnorm(30)), model = c("ANN", "AAN"), control = list(prior_models = c(1, 0)))
  p <- fit$fit$results[[2]]$particles
  expect_null(p$final_mean)
  expect_false(is.null(p$theta))
  expect_true(all(forecast_mixture(fit$fit, 3L)$model == 1))
  expect_false(is.null(summary(fit)$parameters$AAN))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(fit, models = "AAN"))
})

test_that("constant series: exact bounds at h = 1", {
  set.seed(1)
  fit <- suppressWarnings(bets(ts(rep(5, 20), frequency = 4)))
  fc <- predict(fit, h = 3, level = c(80, 95), return_mixture = TRUE)
  df <- fit$control$nu0 + 20
  expect_equal(as.numeric(fc$upper[1, ]), 5 + sqrt(fit$psi0 / df) * stats::qt(c(0.9, 0.975), df),
               tolerance = 1e-10)
  expect_equal(as.numeric(fc$mean), rep(5, 3))
  expect_true(all(fc$mixture$location == 5))
})

test_that("predict.bets() skips failed models", {
  set.seed(1)
  fit <- suppressWarnings(bets(ts(cumsum(rnorm(30))), model = c("ANN", "AAdN"),
                               control = list(min_ess = 1e6, N_iter_max = 1)))
  expect_equal(fit$fit$model_weights[2], 0)
  fc <- predict(fit, h = 4, return_mixture = TRUE)
  expect_true(all(is.finite(c(fc$mean, fc$lower, fc$upper))))
  expect_length(fc$mixture$weights, length(heavy_set(fit$fit$results[[1]]$particles$w, 1e-4)))
})

test_that("t_mixture_quantiles() solves the mixture CDF, also for multimodal mixtures", {
  set.seed(1)
  probs <- c(0.001, 0.025, 0.3, 0.5, 0.9, 0.999)
  for (K in c(1, 3, 50, 2000)) {
    for (df in c(4.5, 30, 300)) {
      w <- stats::rexp(K)
      loc <- matrix(stats::rnorm(K * 3, sd = 5), K)
      sc <- matrix(stats::rexp(K * 3) + 0.1, K)
      q <- t_mixture_quantiles(w, loc, sc, df, probs)
      expect_equal(dim(q), c(3L, length(probs)))
      for (j in 1:3) {
        cdf <- vapply(q[j, ], function(x) sum(w * stats::pt((x - loc[, j]) / sc[, j], df)), numeric(1))
        expect_equal(cdf / sum(w), probs, tolerance = 1e-7, info = sprintf("K = %d, df = %g", K, df))
      }
    }
  }
  # One component: quantiles of the t
  expect_equal(t_mixture_quantiles(1, matrix(3), matrix(2), 7, probs)[1, ], 3 + 2 * stats::qt(probs, 7))
  # Two modes 100 scales apart, quantiles in the gap between them
  loc <- matrix(c(0, 100))
  sc <- matrix(c(1, 2))
  for (w in list(c(0.5, 0.5), c(0.99, 0.01))) {
    for (df in c(4.5, 30)) {
      p <- w[1] + c(-0.005, -1e-4, 0, 1e-4, 0.005)
      q <- t_mixture_quantiles(w, loc, sc, df, p)[1, ]
      cdf <- vapply(q, function(x) sum(w * stats::pt((x - loc) / sc, df)), numeric(1))
      expect_equal(cdf, p, tolerance = 1e-7)
      root <- stats::uniroot(function(x) sum(w * stats::pt((x - loc) / sc, df)) - p[5], c(-50, 150),
                             tol = 1e-12)$root
      expect_equal(q[5], root, tolerance = 1e-6)
    }
  }
  expect_equal(t_mixture_quantiles(1, matrix(0), matrix(1), 5, c(0, 1))[1, ], c(-Inf, Inf))
  expect_error(t_mixture_quantiles(1, matrix(0), matrix(0), 5, 0.5), "positive scales")
})

test_that("closed-form forecasts agree with 1e5 paths simulated from the same mixture", {
  # All six models; horizons beyond two seasonal periods. z: differences in units of
  # their Monte Carlo error (mean; CDF of the paths at the closed-form quantiles)
  skip_on_cran()
  n <- 1e5
  probs <- c(0.025, 0.1, 0.5, 0.9, 0.975)
  all_models <- c("ANN", "AAN", "AAdN", "ANA", "AAA", "AAdA")
  cases <- list(list(y = USAccDeaths, h = 30L, models = all_models),
                list(y = UKgas, h = 10L, models = all_models),
                list(y = Nile, h = 10L, models = c("ANN", "AAN", "AAdN")))
  for (cs in cases) {
    z <- NULL
    for (code in cs$models) {
      set.seed(1)
      fc <- predict(bets(cs$y, model = code), h = cs$h, n_traj = n, return_paths = TRUE,
                    return_mixture = TRUE)
      mix <- fc$mixture
      paths <- unclass(fc$paths)
      q <- t_mixture_quantiles(mix$weights, mix$location, mix$scale, mix$df, probs)
      F_paths <- vapply(seq_len(cs$h), function(j) colMeans(outer(paths[j, ], q[j, ], "<=")),
                        numeric(length(probs)))
      z <- c(z, (F_paths - probs) / sqrt(probs * (1 - probs) / n),
             (rowMeans(paths) - fc$mean) / (apply(paths, 1, stats::sd) / sqrt(n)))
    }
    label <- sprintf("frequency %d", stats::frequency(cs$y))
    expect_lt(max(abs(z)), 5, label = paste("max|z|,", label))
    expect_lt(sqrt(mean(z^2)), 1.5, label = paste("rms(z),", label))
  }
})
