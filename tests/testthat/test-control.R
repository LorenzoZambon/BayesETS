# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() fills in all defaults for empty list", {
  ctrl    <- resolve_bets_control(list())
  defaults <- bets_control_defaults()
  expect_named(ctrl, names(defaults), ignore.order = TRUE)
})

test_that("resolve_bets_control() uses correct default N_iter_max", {
  ctrl <- resolve_bets_control(list())
  expect_equal(ctrl$N_iter_max, 30)
})

test_that("resolve_bets_control() uses correct default lr", {
  ctrl <- resolve_bets_control(list())
  expect_equal(ctrl$lr, 0.9)
})

test_that("resolve_bets_control() uses correct default nu0", {
  ctrl <- resolve_bets_control(list())
  expect_equal(ctrl$nu0, 3)
})

test_that("resolve_bets_control() uses correct default is_df", {
  ctrl <- resolve_bets_control(list())
  expect_equal(ctrl$is_df, 5)
})

# ---------------------------------------------------------------------------
# Overrides
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() overrides N_iter_max", {
  ctrl <- resolve_bets_control(list(N_iter_max = 5))
  expect_equal(ctrl$N_iter_max, 5)
})

test_that("resolve_bets_control() overrides lr while keeping other defaults", {
  ctrl <- resolve_bets_control(list(lr = 0.5))
  expect_equal(ctrl$lr, 0.5)
  expect_equal(ctrl$N_iter_max, 30)
})

test_that("resolve_bets_control() overrides N_draw", {
  ctrl <- resolve_bets_control(list(N_draw = 64L))
  expect_equal(ctrl$N_draw, 64L)
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() rejects non-list input and unnamed entries", {
  expect_error(resolve_bets_control("not_a_list"), "control must be a named list")
  expect_error(resolve_bets_control(list(64L)), "control must be a named list")
})

test_that("resolve_bets_control() rejects nu0 <= 2 (no prior mean of sigma^2)", {
  expect_error(resolve_bets_control(list(nu0 = 2)), "nu0 must be")
  expect_no_error(resolve_bets_control(list(nu0 = 2.5)))
})

test_that("resolve_bets_control() rejects invalid values", {
  bad <- list(N_iter_max = 0, N_draw = c(64, 128.5), N_final = -1, n_quad = 1:5,
              n_scan = -1, min_ess = 0, is_df = 0, is_scale = NA, c_inflate_eta = "3",
              lr = 1.5, nu0 = Inf, psi0 = -1, phi_min = 0, phi_max = 1.1,
              prior_models = c(-1, 2), integration = "mcmc")
  for (nm in names(bad)) {
    expect_error(resolve_bets_control(bad[nm]), paste0("control\\$", nm, " must be"),
                 info = nm)
  }
  expect_error(resolve_bets_control(list(phi_min = 0.9, phi_max = 0.85)),
               "phi_max must be")
})

test_that("resolve_bets_control() accepts valid non-default values", {
  expect_no_error(resolve_bets_control(list(
    N_draw = 512, N_final = c(100, 200), min_ess = c(32, 64), n_scan = 0, psi0 = 2,
    phi_max = 1, lr = 1, prior_models = c(1, 0, 2), integration = "ais")))
})

test_that("resolve_bets_control() rejects unknown keys and lists the valid ones", {
  expect_error(resolve_bets_control(list(foo = 1)), "Unknown control entry 'foo'\\.\nValid entries: N_iter_max")
  expect_error(resolve_bets_control(list(bad_key = 1)), "bad_key")
})

test_that("resolve_bets_control() suggests the closest valid name for a typo", {
  expect_error(resolve_bets_control(list(N_drw = 512)),
               "Unknown control entry 'N_drw'\\. Did you mean 'N_draw'\\?\nValid entries:")
  # Retired names are unknown entries too ('combination' is now a bets() argument)
  expect_error(resolve_bets_control(list(method = "bma")), "Unknown control entry 'method'")
})

test_that("resolve_bets_control() reports every unknown key", {
  expect_error(resolve_bets_control(list(foo = 1, n_qad = 9)),
               "'foo'\\.\nUnknown control entry 'n_qad'\\. Did you mean 'n_quad'\\?")
})

