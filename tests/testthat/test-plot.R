# Shared fitted object, created once for the whole file
local({
  set.seed(1)
  fit_plot <<- bets(ts(10 + rep(c(2, -1, 0, -1), 8) + rnorm(32, sd = 0.5), frequency = 4),
                    model = c("ANN", "ANA", "AAdA"))
})

# Integral of each plotted density over its range
integrals <- function(d) {
  d <- d[!is.na(d$density), ]
  vapply(split(d, list(d$model, d$parameter), drop = TRUE), function(g) {
    sum(diff(g$x) * (g$density[-1] + g$density[-nrow(g)]) / 2)
  }, numeric(1))
}

test_that("plot.bets() draws the shown models and returns their densities", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  d <- plot(fit_plot, ask = FALSE)
  labels <- c("ANN", "ANA", "AAdA")
  expect_setequal(unique(d$model), labels[shown_models(fit_plot$fit$model_weights)])
  expect_named(d, c("model", "parameter", "x", "density", "in_interval"))
  expect_true(all(abs(integrals(d) - 1) < 0.01))
  smoothing <- d$parameter != "sigma"
  expect_true(all(d$x[smoothing] >= 0 & d$x[smoothing] <= 1))
  expect_true(all(d$x[d$parameter == "phi"] >= fit_plot$control$phi_min &
                  d$x[d$parameter == "phi"] <= fit_plot$control$phi_max))
})

test_that("plot.bets(states = TRUE) adds the initial states", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  d <- plot(fit_plot, models = "ANA", states = TRUE, ask = FALSE)
  expect_setequal(unique(d$parameter), c("alpha", "gamma", "sigma", "l", paste0("s", 1:4)))
  expect_true(all(abs(integrals(d) - 1) < 0.01))
})

test_that("plot.bets() plots any fitted model on request, and errors for the others", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  lowest <- c("ANN", "ANA", "AAdA")[which.min(fit_plot$fit$model_weights)]
  expect_no_error(plot(fit_plot, models = lowest, ask = FALSE))
  expect_error(plot(fit_plot, models = "AAA"), "models must be labels of the fitted models")
  set.seed(1)
  f <- suppressWarnings(bets(USAccDeaths, model = c("ANN", "AAdA"),
                             control = list(min_ess = 1e6, N_iter_max = 2)))
  expect_error(plot(f, models = "AAdA"), "No posterior for the failed model")
  expect_error(plot(fit_plot, level = 100), "level must be")
})

test_that("plot.bets() of a constant series shows the level as a point mass", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(1)
  f <- suppressWarnings(bets(ts(rep(5, 20), frequency = 4)))
  d <- plot(f, states = TRUE, ask = FALSE)
  expect_true(is.na(d$density[d$parameter == "l"]) && d$x[d$parameter == "l"] == 5)
})

test_that("season_order(): seasonal states ordered and labelled by season", {
  # Monthly series starting in March: s12 applies to March, s2 to January
  s <- season_order(ts(1:30, start = c(2000, 3), frequency = 12), 12)
  expect_identical(s$labels, month.abb)
  expect_identical(s$states[1:3], c("s2", "s1", "s12"))
  q <- season_order(ts(1:10, start = c(2000, 1), frequency = 4), 4)
  expect_identical(q$labels, paste0("Q", 1:4))
  expect_identical(q$states, c("s4", "s3", "s2", "s1"))
})
