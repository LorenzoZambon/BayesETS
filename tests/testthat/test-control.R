# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() fills in all defaults for empty list", {
  ctrl    <- BETS:::resolve_bets_control(list())
  defaults <- BETS:::bets_control_defaults()
  expect_named(ctrl, names(defaults), ignore.order = TRUE)
})

test_that("resolve_bets_control() uses correct default N_iter_max", {
  ctrl <- BETS:::resolve_bets_control(list())
  expect_equal(ctrl$N_iter_max, 30)
})

test_that("resolve_bets_control() uses correct default lr", {
  ctrl <- BETS:::resolve_bets_control(list())
  expect_equal(ctrl$lr, 0.9)
})

test_that("resolve_bets_control() uses correct default nu0", {
  ctrl <- BETS:::resolve_bets_control(list())
  expect_equal(ctrl$nu0, 3)
})

test_that("resolve_bets_control() uses correct default is_df", {
  ctrl <- BETS:::resolve_bets_control(list())
  expect_equal(ctrl$is_df, 5)
})

# ---------------------------------------------------------------------------
# Overrides
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() overrides N_iter_max", {
  ctrl <- BETS:::resolve_bets_control(list(N_iter_max = 5))
  expect_equal(ctrl$N_iter_max, 5)
})

test_that("resolve_bets_control() overrides lr while keeping other defaults", {
  ctrl <- BETS:::resolve_bets_control(list(lr = 0.5))
  expect_equal(ctrl$lr, 0.5)
  expect_equal(ctrl$N_iter_max, 30)
})

test_that("resolve_bets_control() overrides N_draw", {
  ctrl <- BETS:::resolve_bets_control(list(N_draw = 64L))
  expect_equal(ctrl$N_draw, 64L)
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() rejects non-list input and unnamed entries", {
  expect_error(BETS:::resolve_bets_control("not_a_list"), "control must be a named list")
  expect_error(BETS:::resolve_bets_control(list(64L)), "control must be a named list")
})

test_that("resolve_bets_control() rejects nu0 <= 2 (no prior mean of sigma^2)", {
  expect_error(BETS:::resolve_bets_control(list(nu0 = 2)), "nu0 must be")
  expect_no_error(BETS:::resolve_bets_control(list(nu0 = 2.5)))
})

test_that("resolve_bets_control() rejects unknown keys and lists the valid ones", {
  expect_error(BETS:::resolve_bets_control(list(foo = 1)), "Unknown control entry 'foo'\\.\nValid entries: N_iter_max")
  expect_error(BETS:::resolve_bets_control(list(bad_key = 1)), "bad_key")
})

test_that("resolve_bets_control() suggests the closest valid name for a typo", {
  expect_error(BETS:::resolve_bets_control(list(N_drw = 512)),
               "Unknown control entry 'N_drw'\\. Did you mean 'N_draw'\\?\nValid entries:")
  # Retired names are unknown entries too ('combination' is now a bets() argument)
  expect_error(BETS:::resolve_bets_control(list(method = "bma")), "Unknown control entry 'method'")
})

test_that("resolve_bets_control() reports every unknown key", {
  expect_error(BETS:::resolve_bets_control(list(foo = 1, n_qad = 9)),
               "'foo'\\.\nUnknown control entry 'n_qad'\\. Did you mean 'n_quad'\\?")
})

