#' Predict from a BETS model
#'
#' Forecast distribution of a fitted `bets` object, from simulated
#' trajectories.
#'
#' @param object A fitted object from [bets()].
#' @param newdata Not supported: must be `NULL`.
#' @param h Forecast horizon. Default is the frequency of the series
#'   for seasonal series (if up to 24), and 10 otherwise.
#' @param level Levels of the prediction intervals.
#' @param n_traj Number of simulated trajectories (default 1000).
#'   Increase it to reduce the Monte Carlo error of the interval bounds.
#' @param ... Unused.
#'
#' @return An object of class `bets_forecast`.
#'
#' @examples
#' set.seed(1)
#' fit <- bets(USAccDeaths)
#' fc <- predict(fit)
#' fc
#'
#' # Set forecast horizon and levels of the prediction intervals
#' predict(fit, h = 6, level = c(50, 90))
#'
#' @export
predict.bets <- function(object, newdata = NULL, h = NULL, level = c(80, 95),
                         n_traj = 1000, ...) {
  if (!is.null(newdata))
    stop("`newdata` is not supported for BETS models: forecasts are always ",
         "generated forward from the end of the training series.")

  if (is.null(h)) {                     # default: set h to either 10 or frequency of the series
    m <- round(stats::frequency(object$y))
    h <- if (m > 1 && m <= 24) m else 10
  } else {
    if (!is_count(h))
      stop("h must be a single positive integer")
  }

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

#' @export
print.bets_forecast <- function(x, ...) {
  out <- data.frame(`Point Forecast` = as.numeric(x$mean), check.names = FALSE)
  for (j in seq_along(x$level)) {
    out[[paste("Lo", x$level[j])]] <- as.numeric(x$lower[, j])
    out[[paste("Hi", x$level[j])]] <- as.numeric(x$upper[, j])
  }
  rownames(out) <- time_labels(x$mean)
  print(out, ...)
  invisible(x)
}

# Labels of the time points of a ts: "Jan 2020" (monthly), "2020 Q1" (quarterly),
# the year (yearly), otherwise the time
time_labels <- function(x) {
  t <- as.numeric(stats::time(x))
  f <- stats::frequency(x)
  year <- floor(t + 1e-8)
  period <- round((t - year) * f) + 1
  if (f == 12) {
    paste(month.abb[period], year)
  } else if (f == 4) {
    paste0(year, " Q", period)
  } else if (f == 1) {
    as.character(year)
  } else {
    format(round(t, 3))
  }
}

# Future trajectories (n_traj x h) of the model combination, from the fit
# element of a bets object: each model contributes a number of trajectories
# proportional to its weight
simulate_future_trajectories <- function(bets_fit, h = 10, n_traj = 1000) {
  n_models <- length(bets_fit$results)
  model_weights <- bets_fit$model_weights

  traj_list <- vector("list", n_models)
  n_traj_list <- draws_per_model(model_weights, n_traj)

  for (i in seq_along(bets_fit$results)) {
    if (n_traj_list[i] > 0) {
      idxs <- sample(nrow(bets_fit$results[[i]]$thetas), size = n_traj_list[i], replace = TRUE)
      res_i <- bets_fit$results[[i]]
      traj_list[[i]] <- ets_future_traj(
        model_components = res_i$model_components,
        states = res_i$states[idxs, , drop = FALSE],
        params = res_i$thetas[idxs, , drop = FALSE],
        sigma2s = res_i$sigma2s[idxs],
        h = h
      )
    }
  }

  do.call(rbind, traj_list)
}

# Simulated future trajectories of an ETS model, one per row of params (with
# the final states and \sigma^2 of the same posterior draw)
ets_future_traj <- function(model_components, states, params, sigma2s, h = 10) {
  flags <- model_flags(model_components)
  trend <- flags$trend
  seas <- flags$seas
  damped <- flags$damped

  N_samples <- nrow(params)

  # Seasonal period from the names of the states
  s_cols <- grep("^s\\d+$", colnames(states), value = TRUE)
  if (length(s_cols) > 1) {
    s_idx <- as.integer(sub("^s", "", s_cols))
    s_cols <- s_cols[order(s_idx)]
  }
  m <- if (seas) length(s_cols) else 1
  if (seas && m == 0) stop("Seasonal model but no seasonal states found.")

  alpha <- params[, "alpha"]
  if (trend) {
    beta <- params[, "beta"]
    phi  <- if (damped) params[, "phi"] else 1
  }
  if (seas) gamma <- params[, "gamma"]

  l <- states[, "l"]
  if (trend) b <- states[, "b"]
  if (seas)  s <- states[, s_cols, drop = FALSE]

  # Errors for all horizons (one variance per draw)
  errors <- matrix(stats::rnorm(N_samples * h), nrow = N_samples, ncol = h) * sqrt(sigma2s)

  forecasts <- matrix(nrow = N_samples, ncol = h)
  for (i in 1:h) {
    forecasts[, i] <- l
    if (trend) forecasts[, i] <- forecasts[, i] + phi * b
    if (seas)  forecasts[, i] <- forecasts[, i] + s[, ((i - 1) %% m) + 1]
    forecasts[, i] <- forecasts[, i] + errors[, i]

    # State update
    l <- l + alpha * errors[, i]
    if (trend) {
      l <- l + phi * b
      b <- phi * b + beta * errors[, i]
    }
    if (seas) {
      s[, ((i - 1) %% m) + 1] <- s[, ((i - 1) %% m) + 1] + gamma * errors[, i]
    }
  }

  forecasts
}
