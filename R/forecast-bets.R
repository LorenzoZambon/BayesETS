#' Predict from a BETS model
#'
#' Forecast distribution of a fitted `bets` object, from simulated
#' trajectories.
#'
#' @param object A fitted object from [bets()].
#' @param newdata Not supported: must be `NULL`.
#' @param h Forecast horizon.
#' @param level Levels of the prediction intervals.
#' @param n_traj Number of simulated trajectories (default 1000). With the
#'   default, the Monte Carlo error of the interval bounds is about 2% of the
#'   width of the 95% interval; it decreases as `1 / sqrt(n_traj)`.
#' @param ... Unused.
#'
#' @return An object of class `bets_forecast`.
#' @export
predict.bets <- function(object, newdata = NULL, h = 10, level = c(80, 95),
                         n_traj = 1000, ...) {
  if (!is.null(newdata))
    stop("`newdata` is not supported for BETS models: forecasts are always ",
         "generated forward from the end of the training series.")

  if (!is_count(h))
    stop("h must be a single positive integer")
  if (!is_count(n_traj))
    stop("n_traj must be a single positive integer")

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
  deltat  <- 1 / tspx[3]          # time step
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
