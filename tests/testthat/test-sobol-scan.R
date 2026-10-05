# ---------------------------------------------------------------------------
# sobol_scan_rb(): importance weights of the Sobol early exit
# ---------------------------------------------------------------------------

# ETS(A,N,A) series: for d >= 2 the Sobol proposal is not uniform on the
# admissible region (gamma = (1 - alpha) * u has density 1 / (1 - alpha)).
sim_ana <- function(L, m, alpha, gamma, sigma) {
  l <- 100
  s <- 6 * sin(2 * pi * seq_len(m) / m)
  y <- numeric(L)
  for (t in seq_len(L)) {
    j <- (t - 1) %% m + 1
    e <- stats::rnorm(1, 0, sigma)
    y[t] <- l + s[j] + e
    l <- l + alpha * e
    s[j] <- s[j] + gamma * e
  }
  stats::ts(y, frequency = m)
}

test_that("Sobol early-exit log evidence matches quadrature for d = 2", {
  set.seed(4)
  y <- sim_ana(48, 4, alpha = 0.3, gamma = 0.2, sigma = 1.5)
  mc <- BETS:::coerce_model_components("ANA", 4)
  fit_single <- function(control) {
    set.seed(1)
    ctrl <- BETS:::resolve_bets_control(control, 4)
    BETS:::fit_bets_models(y, mc, ctrl)$results[[1]]
  }
  # min_ess = 1 forces the early exit: the estimate is the Sobol scan alone.
  res_sobol <- fit_single(list(n_sobol = 4096L, min_ess = 1, N_final = 50L))
  res_quad  <- fit_single(list(integration_method = "quadrature", n_quad = 25L, N_final = 50L))
  expect_equal(res_sobol$n_iter, 0L)
  # The weights that treated the Sobol points as prior draws were off by +0.38 here.
  expect_lt(abs(res_sobol$log_evidence - res_quad$log_evidence), 0.02)
})
