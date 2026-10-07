# ---------------------------------------------------------------------------
# 3-character string codes
# ---------------------------------------------------------------------------

test_that("'ANN' parses to correct components", {
  mc <- BETS:::coerce_model_components("ANN", m = 1)
  expect_identical(mc[[1]], c("A", "N", "N", "FALSE"))
})

test_that("'AAN' parses correctly (trend, no season)", {
  mc <- BETS:::coerce_model_components("AAN", m = 1)
  expect_identical(mc[[1]], c("A", "A", "N", "FALSE"))
})

test_that("'AAA' parses correctly (trend + season)", {
  mc <- BETS:::coerce_model_components("AAA", m = 12)
  expect_identical(mc[[1]], c("A", "A", "A", "FALSE"))
})

test_that("'ANA' parses correctly (season only, no trend)", {
  mc <- BETS:::coerce_model_components("ANA", m = 4)
  expect_identical(mc[[1]], c("A", "N", "A", "FALSE"))
})

# ---------------------------------------------------------------------------
# 4-character damped codes
# ---------------------------------------------------------------------------

test_that("'AAdN' parses to damped trend", {
  mc <- BETS:::coerce_model_components("AAdN", m = 1)
  expect_identical(mc[[1]], c("A", "A", "N", "TRUE"))
})

test_that("'AAdA' parses to damped trend with season", {
  mc <- BETS:::coerce_model_components("AAdA", m = 12)
  expect_identical(mc[[1]], c("A", "A", "A", "TRUE"))
})

# ---------------------------------------------------------------------------
# ZZZ model space
# ---------------------------------------------------------------------------

test_that("'ZZZ' with m=1 returns 3 non-seasonal models", {
  mc <- BETS:::coerce_model_components("ZZZ", m = 1)
  expect_length(mc, 3)
  seasons <- vapply(mc, `[[`, character(1), 3)
  expect_true(all(seasons == "N"))
})

test_that("'ZZZ' with m=12 returns 6 models", {
  mc <- BETS:::coerce_model_components("ZZZ", m = 12)
  expect_length(mc, 6)
})

test_that("'ZZZ' leaves out seasonal models when n <= m or m > 24", {
  expect_warning(mc <- BETS:::coerce_model_components("ZZZ", m = 12, n = 12), "not sufficient")
  expect_length(mc, 3)
  expect_warning(mc <- BETS:::coerce_model_components("ZZZ", m = 48, n = 200), "> 24")
  expect_length(mc, 3)
})

labels_of <- function(model, m = 12, n = Inf) {
  vapply(BETS:::coerce_model_components(model, m, n = n), BETS:::ets_label, character(1))
}

test_that("'ZZZ' gives the 6 additive models in a fixed order", {
  expect_identical(labels_of("ZZZ"), c("ANN", "AAN", "AAdN", "ANA", "AAA", "AAdA"))
})

test_that("partial 'Z' codes give all the options of their components", {
  expect_identical(labels_of("AZN"), c("ANN", "AAN", "AAdN"))
  expect_identical(labels_of("ANZ"), c("ANN", "ANA"))
  expect_identical(labels_of("ANZ", m = 1), "ANN")    # no warning for non-seasonal data
  expect_warning(lab <- labels_of("ANZ", n = 10), "not sufficient")
  expect_identical(lab, "ANN")
  expect_identical(labels_of(list("ZZZ")), labels_of("ZZZ"))
})

test_that("a vector or list of codes gives one model per code", {
  expect_identical(labels_of(c("ANN", "AAdN", "ANA")), c("ANN", "AAdN", "ANA"))
  expect_identical(labels_of(list("ANN", "AAdN", "ANA")), c("ANN", "AAdN", "ANA"))
})

test_that("duplicate models are removed with a warning", {
  expect_warning(mc <- BETS:::coerce_model_components(c("ANN", "ann", "A N N"), m = 1),
                 "Duplicate models removed: ANN")
  expect_length(mc, 1)
  expect_warning(lab <- labels_of(list("AZN", "AAdN")), "Duplicate models removed: AAdN")
  expect_identical(lab, c("ANN", "AAN", "AAdN"))
})

# ---------------------------------------------------------------------------
# Seasonal period (rules of forecast::ets)
# ---------------------------------------------------------------------------

test_that("seasonal_period() follows forecast::ets()", {
  expect_identical(BETS:::seasonal_period(ts(1:10, frequency = 12)), 12L)
  expect_identical(BETS:::seasonal_period(ts(1:10, frequency = 1)), 1L)
  expect_warning(m <- BETS:::seasonal_period(ts(1:10, frequency = 0.1)), "below 1")
  expect_identical(m, 1L)
  expect_warning(m <- BETS:::seasonal_period(ts(1:10, frequency = 365.25 / 7)),
                 "Non-integer seasonal period")
  expect_identical(m, 1L)
})

test_that("list of models parses each element independently", {
  mc <- BETS:::coerce_model_components(list("ANN", "AAN"), m = 1)
  expect_length(mc, 2)
  expect_identical(mc[[1]], c("A", "N", "N", "FALSE"))
  expect_identical(mc[[2]], c("A", "A", "N", "FALSE"))
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("seasonal model with m=1 raises an error", {
  expect_error(
    BETS:::coerce_model_components("AAA", m = 1),
    "Seasonal models require frequency"
  )
})

test_that("'MNN' with additive.only=TRUE raises an error", {
  expect_error(
    BETS:::coerce_model_components("MNN", m = 1, additive.only = TRUE),
    "Multiplicative components are not allowed"
  )
})

test_that("invalid trend component raises an error", {
  expect_error(
    BETS:::coerce_model_components("AXN", m = 1),
    "Invalid trend component"
  )
})

test_that("invalid error component raises an error", {
  expect_error(
    BETS:::coerce_model_components("XNN", m = 1),
    "Invalid error component"
  )
})

test_that("invalid season component raises an error", {
  expect_error(
    BETS:::coerce_model_components("ANX", m = 1),
    "Invalid seasonal component"
  )
})

test_that("damping is only allowed for an additive or multiplicative trend", {
  expect_error(BETS:::coerce_model_components("ANdN", m = 1), "Invalid trend component")
  expect_error(BETS:::coerce_model_components("AZdA", m = 12), "Invalid trend component")
})

test_that("codes of invalid length and component vectors raise an error", {
  expect_error(BETS:::coerce_model_components("AN", m = 1), "must have 3 characters")
  expect_error(BETS:::coerce_model_components(c("A", "A", "N"), m = 1), "must have 3 characters")
  expect_error(BETS:::coerce_model_components(1, m = 1), "model must be an ETS code")
})
