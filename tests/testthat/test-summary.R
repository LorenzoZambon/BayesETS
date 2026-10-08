# Shared fitted objects, created once for the whole file
local({
  set.seed(1)
  fit_seas <<- bets(ts(10 + rep(c(2, -1, 0, -1), 8) + rnorm(32, sd = 0.5), frequency = 4),
                    model = c("ANN", "ANA", "AAdA"))
})

# ---------------------------------------------------------------------------
# Weighted particles stored in the fit
# ---------------------------------------------------------------------------

test_that("each model stores its weighted particles", {
  n <- length(fit_seas$y)
  for (r in fit_seas$fit$results) {
    p <- r$particles
    expect_identical(colnames(p$theta), colnames(r$thetas))
    expect_length(p$w, nrow(p$theta))
    expect_length(p$sigma2_scale, nrow(p$theta))
    expect_true(all(p$w > 0))
    expect_equal(sum(p$w), 1, tolerance = 1e-12)
    expect_equal(p$sigma2_df, fit_seas$control$nu0 + n)
  }
})

test_that("posterior means from the particles agree with the posterior draws", {
  for (r in fit_seas$fit$results) {
    s <- posterior_summary(r$particles)
    d <- cbind(r$thetas, sigma = sqrt(r$sigma2s))
    se <- apply(d, 2, stats::sd) / sqrt(nrow(d))
    expect_true(all(abs(s[colnames(d), "mean"] - colMeans(d)) < 5 * se + 1e-6),
                info = ets_label(r$model_components))
  }
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

test_that("weighted_quantile() with equal weights is quantile(type = 5)", {
  set.seed(1)
  x <- rnorm(57)
  probs <- c(0.025, 0.3, 0.5, 0.975)
  expect_equal(weighted_quantile(x, rep(2, 57), probs), unname(stats::quantile(x, probs, type = 5)))
  # zero weights are ignored
  expect_equal(weighted_quantile(c(x, 100), c(rep(1, 57), 0), probs),
               weighted_quantile(x, rep(1, 57), probs))
})

test_that("sigma_summary() is exact for a single scaled inverse chi-squared", {
  scale <- 7; df <- 12; probs <- c(0.025, 0.975)
  s <- sigma_summary(rep(scale, 3), df, rep(1 / 3, 3), probs)
  expect_equal(s[3:4], sqrt(scale / stats::qchisq(1 - probs, df)))
  expect_equal(s[1]^2 + s[2]^2, scale / (df - 2))                   # E[sigma^2]
  set.seed(1)
  expect_equal(s[1], mean(sqrt(scale / stats::rchisq(1e6, df))), tolerance = 1e-3)
})

test_that("sigma_summary() quantiles of a mixture solve the mixture CDF", {
  scale <- c(1, 4, 9); w <- c(0.2, 0.5, 0.3); df <- 10
  q <- sigma_summary(scale, df, w, c(0.1, 0.9))[3:4]
  cdf <- function(s) sum(w * stats::pchisq(scale / s^2, df, lower.tail = FALSE))
  expect_equal(c(cdf(q[1]), cdf(q[2])), c(0.1, 0.9), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# summary.bets()
# ---------------------------------------------------------------------------

test_that("summary.bets() returns the summaries of every model", {
  s <- summary(fit_seas)
  expect_s3_class(s, "summary.bets")
  expect_identical(s$models$model, c("ANN", "ANA", "AAdA"))
  expect_equal(s$models$weight, fit_seas$fit$model_weights)
  expect_identical(rownames(s$parameters$AAdA), c("alpha", "beta", "phi", "gamma", "sigma"))
  expect_identical(colnames(s$parameters$ANA), c("mean", "sd", "2.5%", "97.5%"))
  for (p in s$parameters) {
    expect_true(all(p$sd >= 0))
    expect_true(all(p[["2.5%"]] <= p$mean & p$mean <= p[["97.5%"]]))
  }
  expect_named(s$probabilities, c("trend", "damped trend", "seasonality"))
  expect_identical(colnames(summary(fit_seas, level = 80)$parameters$ANN), c("mean", "sd", "10%", "90%"))
})

test_that("summary.bets() errors for an invalid level", {
  expect_error(summary(fit_seas, level = 0.95 * 200), "level must be a single number")
  expect_error(summary(fit_seas, level = c(80, 95)), "level must be a single number")
})

test_that("print.summary.bets() shows a table per shown model, in the order of print()", {
  out <- capture.output(print(summary(fit_seas)))
  shown <- shown_models(fit_seas$fit$model_weights)
  labels <- c("ANN", "ANA", "AAdA")[shown]
  headers <- grep("^  [A-Za-z]+: weight ", out, value = TRUE)
  expect_identical(sub(":.*", "", trimws(headers)), labels)
  expect_true(any(grepl("^  +mean +sd +2\\.5% +97\\.5%$", out)))
  expect_true(any(grepl("^  sigma ", out)))
})

test_that("summary of a constant series: uniform prior of alpha, no log evidence", {
  set.seed(1)
  s <- summary(suppressWarnings(bets(ts(rep(5, 20), frequency = 4))), states = TRUE)
  expect_equal(s$parameters$ANN["alpha", "mean"], 0.5, tolerance = 1e-12)
  expect_equal(unlist(s$parameters$ANN["alpha", 3:4]), c(0.025, 0.975), tolerance = 0.01,
               ignore_attr = TRUE)
  expect_equal(unlist(s$states$ANN["l", ]), c(5, 0, 5, 5), ignore_attr = TRUE)
  out <- capture.output(print(s))
  expect_true(any(grepl("ANN: weight 1.000$", out)))
})

# ---------------------------------------------------------------------------
# Initial states
# ---------------------------------------------------------------------------

test_that("the posterior evaluated again at the particles is the one of the fit", {
  for (i in seq_along(fit_seas$fit$results)) {
    r <- fit_seas$fit$results[[i]]
    ctrl <- fit_seas$control
    ctrl$psi0 <- fit_seas$psi0
    y <- stats::ts(as.numeric(fit_seas$y), frequency = fit_seas$period)
    y <- y - centring_shift(y)
    tn <- colnames(r$particles$theta)
    log_g <- make_log_g_rb(y, r$model_components, tn, ctrl, init_rb_prior(y, r$model_components, tn, ctrl))
    z <- transform_theta_to_unconstrained(r$particles$theta, ctrl$phi_min, ctrl$phi_max)
    expect_equal(as.numeric(log_g(z)$ml_res$posterior_scale), r$particles$sigma2_scale,
                 tolerance = 1e-10)
  }
})

test_that("summary(states = TRUE): normalised initial states agree with the posterior draws", {
  s <- summary(fit_seas, states = TRUE)
  expect_identical(rownames(s$states$ANN), "l")
  expect_identical(rownames(s$states$AAdA), c("l", "b", "s1", "s2", "s3", "s4"))
  for (i in seq_along(fit_seas$fit$results)) {
    r <- fit_seas$fit$results[[i]]
    st <- s$states[[i]]
    # draws normalised as in summary(): seasonal states summing to zero
    e <- r$etas
    sc <- grep("^s", colnames(e))
    if (length(sc) > 0) {
      ms <- rowMeans(e[, sc, drop = FALSE])
      e[, "l"] <- e[, "l"] + ms
      e[, sc] <- e[, sc] - ms
      expect_equal(sum(st[sc, "mean"]), 0, tolerance = 1e-8)
    }
    se <- apply(e, 2, stats::sd) / sqrt(nrow(e))
    expect_true(all(abs(st[colnames(e), "mean"] - colMeans(e)) < 5 * se),
                info = ets_label(r$model_components))
  }
  expect_null(summary(fit_seas)$states)
  expect_error(summary(fit_seas, states = NA), "states must be TRUE or FALSE")
  out <- capture.output(print(s))
  expect_identical(sum(out == "  initial states:"), length(shown_models(fit_seas$fit$model_weights)))
})

test_that("t_mixture_summary() is exact for a single t", {
  probs <- c(0.025, 0.975)
  s <- t_mixture_summary(rep(3, 4), rep(2, 4), 7, rep(0.25, 4), probs)
  expect_equal(s, c(3, 2 * sqrt(7 / 5), 3 + 2 * stats::qt(probs, 7)))
})
