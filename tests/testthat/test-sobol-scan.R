# ---------------------------------------------------------------------------
# sobol_scan_rb(): importance weights of the Sobol scan
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

test_that("Sobol scan log evidence matches quadrature for d = 2", {
  set.seed(4)
  y <- sim_ana(48, 4, alpha = 0.3, gamma = 0.2, sigma = 1.5)
  mc <- BETS:::coerce_model_components("ANA", 4)
  ctrl <- BETS:::resolve_bets_control(list(n_quad = 25L, N_final = 50L), 4)
  res_quad <- BETS:::fit_bets_models(y, mc, ctrl, integration = "quadrature")$results[[1]]

  # Same psi0 as fit_bets_models() uses by default.
  ctrl$psi0 <- BETS:::default_psi0(y)
  theta_names <- c("alpha", "gamma")
  prior <- BETS:::init_rb_prior(y, mc[[1]], theta_names, ctrl)
  set.seed(1)
  scan <- BETS:::sobol_scan_rb(y, mc[[1]], theta_names, ctrl$phi_min, ctrl$phi_max,
                               n_sobol = 4096L, eta0 = prior$eta0, V0 = prior$V0,
                               nu0 = ctrl$nu0, psi0 = prior$psi0, L = length(y),
                               log_prior_theta_const = prior$log_prior_theta_const)
  lw_max <- max(scan$log_w)
  log_evidence_scan <- lw_max + log(mean(exp(scan$log_w - lw_max)))
  # The weights that treated the Sobol points as prior draws were off by +0.38 here.
  expect_lt(abs(log_evidence_scan - res_quad$log_evidence), 0.02)
})
