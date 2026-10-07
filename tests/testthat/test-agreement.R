# ---------------------------------------------------------------------------
# Agreement of quadrature and default-budget AIS with a high-budget AIS
# reference.
#
# Slow (high-budget AIS reference fits), so skipped on CRAN.
# ---------------------------------------------------------------------------

# Simulate an additive ETS series with the components of `code`.
sim_ets <- function(code, freq, seed) {
  set.seed(seed)
  mc <- BETS:::coerce_model_components(code, freq)[[1]]
  trend  <- mc[[2]] == "A"
  seas   <- mc[[3]] == "A"
  damped <- mc[[4]] == "TRUE"
  L <- c(`1` = 40, `4` = 48, `12` = 96)[[as.character(freq)]]

  alpha <- stats::runif(1, 0.1, 0.6)
  beta  <- if (trend) stats::runif(1, 0.02, 0.5) * alpha else 0
  gamma <- if (seas) stats::runif(1, 0.05, 0.5) * (1 - alpha) else 0
  phi   <- if (damped) stats::runif(1, 0.85, 0.95) else 1

  l <- 100
  b <- if (trend) 1 else 0
  s <- if (seas) 6 * sin(2 * pi * seq_len(freq) / freq) else rep(0, freq)
  y <- numeric(L)
  for (t in seq_len(L)) {
    j <- (t - 1) %% freq + 1
    e <- stats::rnorm(1, 0, 1.5)
    y[t] <- l + phi * b + s[j] + e
    l <- l + phi * b + alpha * e
    b <- phi * b + beta * e
    s[j] <- s[j] + gamma * e
  }
  stats::ts(y, frequency = freq)
}

# Fit a single model specification through fit_bets_models() (which resolves psi0).
fit_single <- function(y, code, control, integration = "ais") {
  freq <- stats::frequency(y)
  ctrl <- BETS:::resolve_bets_control(control)
  BETS:::fit_bets_models(y, BETS:::coerce_model_components(code, freq), ctrl,
                         integration = integration)$results[[1]]
}

# Five series per frequency, each fitted with its data-generating model;
# together they cover every theta dimension d = 1..4.
agreement_dgp <- list(
  `1`  = c("ANN", "AAN", "AAdN", "AAN", "AAdN"),
  `4`  = c("ANN", "ANA", "AAA", "AAdN", "AAdA"),
  `12` = c("ANN", "ANA", "AAA", "AAdN", "AAdA")
)
agreement_series <- list()
for (freq in c(1, 4, 12)) {
  for (i in 1:5) {
    code <- agreement_dgp[[as.character(freq)]][i]
    agreement_series[[length(agreement_series) + 1]] <- list(
      y = sim_ets(code, freq, seed = 1000 * freq + i),
      code = code,
      seed = i,
      label = sprintf("%s (frequency %d, series %d)", code, freq, i)
    )
  }
}

ais_ref_ctrl <- list(N_draw = 8192L, min_ess = 4096, N_final = 5000L)
default_ctrl <- list(N_final = 5000L)

# The AIS reference fits are expensive: computed once per series, shared by the tests.
ais_ref_cache <- new.env()
ais_reference <- function(s) {
  if (is.null(ais_ref_cache[[s$label]])) {
    set.seed(s$seed)
    ais_ref_cache[[s$label]] <- fit_single(s$y, s$code, ais_ref_ctrl)
  }
  ais_ref_cache[[s$label]]
}

# ---------------------------------------------------------------------------
# Log evidence and posterior mean of alpha
# ---------------------------------------------------------------------------

test_that("quadrature agrees with AIS where it is the default (d <= 2)", {
  # "auto" uses quadrature only for d <= 2 (biased for d >= 3)
  skip_on_cran()
  for (s in agreement_series) {
    mc <- BETS:::coerce_model_components(s$code, stats::frequency(s$y))[[1]]
    if (BETS:::resolve_integration("auto", mc) != "quadrature") next
    res_ais  <- ais_reference(s)
    set.seed(s$seed)
    res_quad <- fit_single(s$y, s$code, default_ctrl, integration = "quadrature")
    expect_lt(abs(res_quad$log_evidence - res_ais$log_evidence), 0.05,
              label = paste("|log evidence difference| for", s$label))
    expect_lt(abs(mean(res_quad$thetas[, "alpha"]) - mean(res_ais$thetas[, "alpha"])), 0.01,
              label = paste("|posterior mean alpha difference| for", s$label))
  }
})

test_that("default AIS agrees with the reference on log evidence and posterior mean of alpha", {
  # AIS is random: median difference over 3 seeds
  skip_on_cran()
  for (s in agreement_series) {
    res_ref <- ais_reference(s)
    runs <- lapply(1:3, function(k) {
      set.seed(100 * k + s$seed)
      fit_single(s$y, s$code, default_ctrl)
    })
    d_le    <- vapply(runs, function(r) r$log_evidence - res_ref$log_evidence, numeric(1))
    d_alpha <- vapply(runs, function(r) mean(r$thetas[, "alpha"]) - mean(res_ref$thetas[, "alpha"]),
                      numeric(1))
    expect_lt(abs(stats::median(d_le)), 0.05,
              label = paste("|median log evidence difference| for", s$label))
    expect_lt(abs(stats::median(d_alpha)), 0.01,
              label = paste("|median posterior mean alpha difference| for", s$label))
  }
})

# ---------------------------------------------------------------------------
# Node-count stability
# ---------------------------------------------------------------------------

test_that("quadrature log evidence is stable between n_quad = 5 and n_quad = 7", {
  skip_on_cran()
  for (s in agreement_series) {
    set.seed(s$seed)
    le5 <- fit_single(s$y, s$code, list(n_quad = 5L, N_final = 50L),
                      integration = "quadrature")$log_evidence
    le7 <- fit_single(s$y, s$code, list(n_quad = 7L, N_final = 50L),
                      integration = "quadrature")$log_evidence
    expect_lt(abs(le7 - le5), 0.05,
              label = paste("|log evidence(n_quad = 7) - log evidence(n_quad = 5)| for", s$label))
  }
})

# ---------------------------------------------------------------------------
# Forecast quantiles
# ---------------------------------------------------------------------------

test_that("80% and 95% forecast quantiles agree with the reference within Monte Carlo error", {
  skip_on_cran()
  n_traj <- 20000L
  probs  <- c(0.025, 0.1, 0.9, 0.975)
  for (freq in c(1, 4, 12)) {
    y <- sim_ets(if (freq == 1) "AAN" else "AAA", freq, seed = 7 * freq)
    h <- max(6, 2 * freq)

    set.seed(1)
    fit_ais  <- bets(y, control = list(integration = "ais", N_final = 4000L,
                                       N_draw = 4096L, min_ess = 2048))
    traj_ais <- BETS:::simulate_future_trajectories(fit_ais$fit, h = h, n_traj = n_traj)

    for (method in c("auto", "quadrature", "ais")) {
      set.seed(1)
      fit <- bets(y, control = list(integration = method, N_final = 4000L))
      fc <- predict(fit, h = h, level = c(80, 95), n_traj = n_traj)
      q_fit <- rbind(as.numeric(fc$lower[, "95"]), as.numeric(fc$lower[, "80"]),
                     as.numeric(fc$upper[, "80"]), as.numeric(fc$upper[, "95"]))

      # Level of each quantile under the AIS predictive, standardised by the
      # binomial error (AIS vs AIS: rms(z) 0.8-1.3, max|z| < 4)
      F_ais <- vapply(seq_len(h), function(j) colMeans(outer(traj_ais[, j], q_fit[, j], "<=")),
                      numeric(length(probs)))
      z <- (F_ais - probs) / sqrt(2 * probs * (1 - probs) / n_traj)
      expect_lt(sqrt(mean(z^2)), 2, label = sprintf("rms(z) for %s, frequency %d", method, freq))
      expect_lt(max(abs(z)), 5, label = sprintf("max|z| for %s, frequency %d", method, freq))
    }
  }
})
