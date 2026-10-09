# Benchmark of the README: bets() vs forecast::ets() (additive models) on all the yearly,
# quarterly and monthly M3 series, with the M3 test sets and horizons.
# Writes benchmark/m3_results.csv, one row per series and method.
# Run from the package root, with BayesETS installed (a few minutes):
#   Rscript benchmark/m3_benchmark.R

library(BayesETS)
library(forecast)
library(Mcomp)

set.seed(123)
periods <- c("YEARLY", "QUARTERLY", "MONTHLY")

# Scale of MASE and MSIS: in-sample MAE of the seasonal naive forecast
scale_of <- function(x) mean(abs(diff(as.numeric(x), lag = frequency(x))))

# Interval score of the (1 - a) intervals, averaged over the horizon
interval_score <- function(y, lo, hi, a = 0.05) {
  mean((hi - lo) + 2 / a * (lo - y) * (y < lo) + 2 / a * (y - hi) * (y > hi))
}

scores <- function(series, period, method, y, fc_mean, lo80, hi80, lo95, hi95, scale, time) {
  data.frame(series = series, period = period, method = method, h = length(y),
             MASE  = mean(abs(y - fc_mean)) / scale,
             cov80 = mean(y >= lo80 & y <= hi80),
             cov95 = mean(y >= lo95 & y <= hi95),
             MSIS  = interval_score(y, lo95, hi95) / scale,
             time  = time)
}

elapsed <- function() proc.time()[["elapsed"]]

rows <- list()
for (p in periods) {
  for (s in Filter(function(s) s$period == p, M3)) {
    h <- length(s$xx)
    y <- as.numeric(s$xx)
    sc <- scale_of(s$x)

    t0 <- elapsed()
    fc <- predict(suppressWarnings(bets(s$x)), h = h, level = c(80, 95))
    t_bets <- elapsed() - t0

    t0 <- elapsed()
    efc <- forecast(ets(s$x, additive.only = TRUE), h = h, level = c(80, 95))
    t_ets <- elapsed() - t0

    rows[[length(rows) + 1]] <- scores(s$sn, p, "bets", y, as.numeric(fc$mean),
                                       fc$lower[, "80"], fc$upper[, "80"],
                                       fc$lower[, "95"], fc$upper[, "95"], sc, t_bets)
    rows[[length(rows) + 1]] <- scores(s$sn, p, "ets", y, as.numeric(efc$mean),
                                       efc$lower[, "80%"], efc$upper[, "80%"],
                                       efc$lower[, "95%"], efc$upper[, "95%"], sc, t_ets)
  }
}

write.csv(do.call(rbind, rows), "benchmark/m3_results.csv", row.names = FALSE)
