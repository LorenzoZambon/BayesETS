#' Forecast from a BETS model
#'
#' Produces forecast distributions from a fitted `bets` object using
#' Monte Carlo trajectories.
#'
#' @param object A fitted object from [bets()].
#' @param h Forecast horizon.
#' @param level Confidence levels for intervals.
#' @param n_traj Number of simulated trajectories.
#' @param ... Unused.
#'
#' @return An object with class `c("bets_forecast", "forecast")`.
#' @export
forecast.bets <- function(object, h = 10, level = c(80, 95), n_traj = NULL, ...) {
  if (is.null(n_traj)) {
    n_traj <- object$control$n_traj_forecast
  }
  if (!is.numeric(h) || length(h) != 1 || h <= 0) {
    stop("h must be a single positive number")
  }

  level <- sort(unique(as.numeric(level)))
  if (any(level <= 0 | level >= 100)) {
    stop("level must be between 0 and 100")
  }

  traj <- simulate_future_trajectories(object$fit, h = as.integer(h), n_traj = as.integer(n_traj))
  mean_fc <- colMeans(traj)

  lower <- sapply(level, function(lv) {
    apply(traj, 2, stats::quantile, probs = (100 - lv) / 200, na.rm = TRUE)
  })
  upper <- sapply(level, function(lv) {
    apply(traj, 2, stats::quantile, probs = 1 - (100 - lv) / 200, na.rm = TRUE)
  })

  if (is.null(dim(lower))) {
    lower <- matrix(lower, ncol = 1)
    upper <- matrix(upper, ncol = 1)
  }
  colnames(lower) <- as.character(level)
  colnames(upper) <- as.character(level)

  x <- object$y
  tspx <- stats::tsp(x)
  deltat <- tspx[3]
  start_fc <- tspx[2] + deltat

  mean_ts <- stats::ts(mean_fc, start = start_fc, frequency = stats::frequency(x))
  lower_ts <- stats::ts(lower, start = start_fc, frequency = stats::frequency(x))
  upper_ts <- stats::ts(upper, start = start_fc, frequency = stats::frequency(x))

  out <- list(
    method = "BETS",
    model = object,
    mean = mean_ts,
    lower = lower_ts,
    upper = upper_ts,
    level = level,
    x = x,
    series = deparse(substitute(object))
  )
  class(out) <- c("bets_forecast", "forecast")
  out
}
