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

test_that("resolve_bets_control() uses correct default eta_df", {
  ctrl <- BETS:::resolve_bets_control(list())
  expect_equal(ctrl$eta_df, 7)
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

test_that("resolve_bets_control() overrides n_sobol", {
  ctrl <- BETS:::resolve_bets_control(list(n_sobol = 128L))
  expect_equal(ctrl$n_sobol, 128L)
})

# ---------------------------------------------------------------------------
# Allowed special keys (method, sampler are stripped, not stored)
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() allows 'method' key without error", {
  expect_no_error(BETS:::resolve_bets_control(list(method = "bma")))
})

test_that("resolve_bets_control() allows 'sampler' key without error", {
  expect_no_error(BETS:::resolve_bets_control(list(sampler = "ais")))
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("resolve_bets_control() rejects non-list input", {
  expect_error(BETS:::resolve_bets_control("not_a_list"), "control must be a named list")
})

test_that("resolve_bets_control() rejects unknown keys", {
  expect_error(BETS:::resolve_bets_control(list(foo = 1)), "Unknown control entries")
})

test_that("resolve_bets_control() reports the unknown key name in the error", {
  expect_error(BETS:::resolve_bets_control(list(bad_key = 1)), "bad_key")
})

test_that("resolve_bets_control() rejects multiple unknown keys", {
  expect_error(BETS:::resolve_bets_control(list(foo = 1, bar = 2)), "Unknown control entries")
})
