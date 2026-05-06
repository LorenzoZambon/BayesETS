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
# 3/4-element character vectors
# ---------------------------------------------------------------------------

test_that("3-element vector c('A','A','N') parses correctly", {
  mc <- BETS:::coerce_model_components(c("A", "A", "N"), m = 1)
  expect_identical(mc[[1]], c("A", "A", "N", "FALSE"))
})

test_that("4-element vector with damped='TRUE' parses correctly", {
  mc <- BETS:::coerce_model_components(c("A", "A", "N", "TRUE"), m = 1)
  expect_identical(mc[[1]], c("A", "A", "N", "TRUE"))
})

test_that("4-element vector with logical damped=TRUE parses correctly", {
  mc <- BETS:::coerce_model_components(c("A", "A", "N", "TRUE"), m = 1)
  expect_identical(mc[[1]], c("A", "A", "N", "TRUE"))
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

test_that("damped trend with no trend component raises an error", {
  expect_error(
    BETS:::coerce_model_components(c("A", "N", "N", "TRUE"), m = 1),
    "Damped trend is only valid"
  )
})

test_that("invalid damped value raises an error", {
  expect_error(
    BETS:::coerce_model_components(c("A", "A", "N", "MAYBE"), m = 1),
    "Invalid damped component"
  )
})
