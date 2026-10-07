#' Predict from a BETS model
#'
#' Produces forecast distributions from a fitted `bets` object using
#' Monte Carlo trajectories.
#'
#' @param object A fitted object from [bets()].
#' @param newdata Not used. Passing anything other than `NULL` raises an error,
#'   since BETS forecasts forward from the end of the training series.
#' @param h Forecast horizon.
#' @param level Confidence levels for intervals.
#' @param n_traj Number of simulated trajectories (default 1000). The Monte
#'   Carlo error of the interval bounds decreases as `1 / sqrt(n_traj)`: with
#'   the default it is about 2% of the width of the 95% interval, and 4 times
#'   more trajectories halve it (the cost grows linearly with `n_traj * h`).
#' @param ... Unused.
#'
#' @return An object of class `bets_forecast`.
#' @export
predict.bets <- function(object, newdata = NULL, h = 10, level = c(80, 95),
                         n_traj = 1000, ...) {
  if (!is.null(newdata))
    stop("`newdata` is not supported for BETS models: forecasts are always ",
         "generated forward from the end of the training series.")

  if (!is.numeric(h) || length(h) != 1 || h <= 0)
    stop("h must be a single positive number")
  if (!is.numeric(n_traj) || length(n_traj) != 1 || n_traj < 1)
    stop("n_traj must be a single positive number")

  level <- sort(unique(as.numeric(level)))
  if (any(level <= 0 | level >= 100))
    stop("level must be between 0 and 100")

  traj <- simulate_future_trajectories(object$fit, h = as.integer(h),
                                       n_traj = as.integer(n_traj))
  mean_fc <- colMeans(traj)

  lower <- sapply(level, function(lv)
    apply(traj, 2, stats::quantile, probs = (100 - lv) / 200, na.rm = TRUE))
  upper <- sapply(level, function(lv)
    apply(traj, 2, stats::quantile, probs = 1 - (100 - lv) / 200, na.rm = TRUE))

  # sapply() drops to a vector when h or the number of levels is 1
  lower <- matrix(lower, nrow = h)
  upper <- matrix(upper, nrow = h)
  colnames(lower) <- as.character(level)
  colnames(upper) <- as.character(level)

  x       <- object$y
  tspx    <- stats::tsp(x)
  deltat  <- 1 / tspx[3]          # tspx[3] is frequency; time step = 1/frequency
  start_fc <- tspx[2] + deltat

  mean_ts  <- stats::ts(mean_fc, start = start_fc, frequency = stats::frequency(x))
  lower_ts <- stats::ts(lower,   start = start_fc, frequency = stats::frequency(x))
  upper_ts <- stats::ts(upper,   start = start_fc, frequency = stats::frequency(x))

  out <- list(
    method = "BETS",
    model  = object,
    mean   = mean_ts,
    lower  = lower_ts,
    upper  = upper_ts,
    level  = level,
    x      = x,
    series = deparse(substitute(object))
  )
  class(out) <- "bets_forecast"
  out
}
