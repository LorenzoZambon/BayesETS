#' Predict from a BETS model
#'
#' Forecast distribution of a fitted `bets` object, computed in closed form.
#' Given the smoothing parameters of a model, the forecast at each horizon has a
#' Student-t distribution, since the initial states and the error variance are
#' integrated analytically. Over the posterior particles of the smoothing
#' parameters (quadrature nodes or AIS draws) and over the models, the forecast
#' distribution is thus a mixture of Student-t distributions, whose quantiles
#' give the prediction intervals. Given the fit, the forecasts are deterministic.
#' Components with negligible weight are left out: the models with weight below
#' 0.0001, and the smallest components holding 0.01% of the weight.
#' Sample paths are simulated only on request.
#'
#' @param object A fitted object from [bets()].
#' @param newdata Not supported: must be `NULL`.
#' @param h Forecast horizon. Default is the frequency of the series
#'   for seasonal series (if up to 24), and 10 otherwise.
#' @param level Levels of the prediction intervals.
#' @param return_paths Logical: also return sample paths from the forecast
#'   distribution (default `FALSE`).
#' @param n_traj Number of sample paths (default 1000), used only if
#'   `return_paths = TRUE`.
#' @param return_mixture Logical: also return the forecast distribution
#'   as a mixture of Student-t distributions, horizon by horizon
#'   (default `FALSE`). The mixture describes each horizon separately:
#'   for quantities involving several horizons (e.g. the total over the
#'   next 12 months), use the sample paths (`return_paths = TRUE`).
#' @param ... Unused.
#'
#' @return An object of class `bets_forecast`, with `print()` and `plot()`
#'   methods ([plot.bets_forecast()]): a list with the point forecasts `mean`
#'   and the interval bounds `lower` and `upper` (one column per level), as time
#'   series after the end of the series, the `level`s, the series `x`, the
#'   fitted `model` and the `method` ("BETS"). With `return_paths = TRUE`, also
#'   `paths`: the sample paths as a time series with `h` rows and `n_traj`
#'   columns (independent draws from the forecast distribution). With
#'   `return_mixture = TRUE`, also `mixture`: a list with the `weights` of the
#'   components (summing to 1), their degrees of freedom `df`, and their
#'   `location` and `scale` (one row per component, one column per horizon).
#'   The forecast CDF at horizon `j` is
#'   `sum(weights * pt((x - location[, j]) / scale[, j], df))`.
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
#' plot(predict(fit, h = 24))
#'
#' # Sample paths, e.g. for the forecast distribution of the total over 12 months
#' fc <- predict(fit, h = 12, return_paths = TRUE)
#' quantile(colSums(fc$paths), c(0.05, 0.5, 0.95))
#'
#' # Probability that the next observation exceeds 8000, from the mixture
#' mix <- predict(fit, h = 1, return_mixture = TRUE)$mixture
#' sum(mix$weights * pt((8000 - mix$location[, 1]) / mix$scale[, 1], mix$df,
#'                      lower.tail = FALSE))
#'
#' @export
predict.bets <- function(object, newdata = NULL, h = NULL, level = c(80, 95),
                         return_paths = FALSE, n_traj = 1000, return_mixture = FALSE, ...) {
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

  if (!(is.logical(return_paths) && length(return_paths) == 1 && !is.na(return_paths)))
    stop("return_paths must be TRUE or FALSE")

  if (!is_count(n_traj))
    stop("n_traj must be a single positive integer")

  if (!(is.logical(return_mixture) && length(return_mixture) == 1 && !is.na(return_mixture)))
    stop("return_mixture must be TRUE or FALSE")

  level <- sort(unique(as.numeric(level)))
  if (any(level <= 0 | level >= 100))
    stop("level must be between 0 and 100")

  h <- as.integer(h)
  mix <- forecast_mixture(object$fit, h)
  mean_fc <- colSums(mix$weights * mix$location)

  # Interval bounds: quantiles of the mixture, all horizons and levels at once (C++)
  a <- (100 - level) / 200
  q <- t_mixture_quantiles(mix$weights, mix$location, mix$scale, mix$df, c(a, 1 - a))
  lower <- q[, seq_along(level), drop = FALSE]
  upper <- q[, length(level) + seq_along(level), drop = FALSE]
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
  if (return_paths) {
    out$paths <- stats::ts(t(simulate_paths(object$fit, mix, h, as.integer(n_traj))),
                           start = start_fc, frequency = stats::frequency(x))
  }
  if (return_mixture) out$mixture <- mix[c("weights", "df", "location", "scale")]
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

#' Plot BETS forecasts
#'
#' Plots the series, the point forecasts and the prediction intervals
#' (the higher the level, the lighter its shade). 
#'
#' @param x A forecast object from [predict.bets()].
#' @param include Number of past observations of the series to plot.
#'   Default: all.
#' @param main Title of the plot. Default: none.
#' @param ... Graphical parameters passed to [graphics::plot()], e.g. `xlab`,
#'   `ylab` or `ylim`.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' set.seed(1)
#' fit <- bets(USAccDeaths)
#' plot(predict(fit, h = 24), include = 36)
#'
#' @export
plot.bets_forecast <- function(x, include = NULL, main = NULL, ...) {
  y <- x$x
  if (!is.null(include)) {
    if (!is_count(include)) stop("include must be a single positive integer")
    if (include < length(y)) y <- stats::window(y, start = stats::time(y)[length(y) - include + 1])
  }
  t_fc <- as.numeric(stats::time(x$mean))
  args <- utils::modifyList(
    list(x = NA, type = "n", main = main, xlab = "", ylab = "",
         xlim = range(stats::time(y), t_fc), ylim = range(y, x$lower, x$upper)),
    list(...))
  do.call(graphics::plot, args)

  # Intervals from the widest to the narrowest, each one covering the previous
  col <- "steelblue"
  graphics::lines(y)
  for (j in rev(seq_along(x$level))) {
    shade <- interval_shade(x$level[j], col)
    if (length(t_fc) > 1) {
      graphics::polygon(c(t_fc, rev(t_fc)), c(x$lower[, j], rev(x$upper[, j])), col = shade, border = NA)
    } else {
      graphics::segments(t_fc, x$lower[, j], t_fc, x$upper[, j], col = shade, lwd = 8, lend = 1)
    }
  }
  if (length(t_fc) > 1) {
    graphics::lines(t_fc, x$mean, col = col, lwd = 2)
  } else {
    graphics::points(t_fc, x$mean, col = col, pch = 19)
  }
  invisible(x)
}

# Opaque shade of a prediction interval, between white and col: lighter for
# higher levels, and the same for the same level in every plot
interval_shade <- function(level, col) {
  ramp <- grDevices::colorRamp(c("white", col))
  grDevices::rgb(ramp(0.15 + 0.75 * sqrt(1 - level / 100)), maxColorValue = 255)
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


################################################################################
# FORECAST DISTRIBUTION

# The smallest components holding together at most .TRUNCATION_MASS of the
# weight are left out of the forecast distribution
.TRUNCATION_MASS <- 1e-4

# Forecast distribution of the model combination, horizon by horizon: mixture
# of t distributions with df degrees of freedom, one per model and particle,
# with weights (model weight x particle weight), location and scale
# (components x h), and the model and particle of each component
forecast_mixture <- function(bets_fit, h) {
  results <- bets_fit$results
  # Models with final states (none for failed and negligible models). Each one
  # contributes its own t components. A model without closed form (e.g.
  # multiplicative ETS) would contribute weighted sample paths instead: the CDF
  # solved in t_mixture_quantiles() would then add the weighted empirical CDF of
  # the paths, and bisection steps would find the quantiles in the same bracket.
  models <- which(vapply(results, function(r) !is.null(r$particles$final_mean), logical(1)))
  n_k <- vapply(results[models], function(r) length(r$particles$w), integer(1))
  model <- rep(models, n_k)
  particle <- sequence(n_k)
  w <- unlist(lapply(models, function(i) bets_fit$model_weights[i] * results[[i]]$particles$w))

  keep <- heavy_set(w, .TRUNCATION_MASS)
  model <- model[keep]
  particle <- particle[keep]
  parts <- lapply(unique(model), function(i) ets_t_mixture(results[[i]], h, particle[model == i]))
  list(weights = w[keep] / sum(w[keep]),
       df = results[[models[1]]]$particles$sigma2_df,
       location = do.call(rbind, lapply(parts, `[[`, "location")),
       scale = do.call(rbind, lapply(parts, `[[`, "scale")),
       model = model, particle = particle)
}

# t components of the forecast distribution of an additive ETS model, at its
# particles idx: given theta, y_{T+j} ~ t_df(f_j' m_T, (sigma2_scale / df) (v_j + f_j' P_T f_j)),
# with m_T = final_mean, P_T = t(final_R) final_R, f_j' x_T = l + \phi_j b + (seasonal
# state of T + j) the forecast without errors and \sigma^2 v_j the variance of the errors
# (C++); list(location, scale), components x h
ets_t_mixture <- function(result, h, idx) {
  p <- result$particles
  flags <- model_flags(result$model_components)
  theta <- p$theta[idx, , drop = FALSE]
  par <- function(name, value) if (name %in% colnames(theta)) theta[, name] else rep(value, length(idx))
  ets_t_components(p$final_mean, p$final_R, idx, theta[, "alpha"], par("beta", 0), par("phi", 1),
                   par("gamma", 0), p$sigma2_scale[idx] / p$sigma2_df, flags$trend, flags$seas, h)
}

# Sample paths (n_traj x h) from the mixture: for each path a component (model
# and particle), \sigma^2 and the final states from their posteriors, then the
# errors and the recursion of the model
simulate_paths <- function(bets_fit, mix, h, n_traj) {
  k <- sample.int(length(mix$weights), n_traj, replace = TRUE, prob = mix$weights)
  paths <- matrix(0, n_traj, h)
  for (i in unique(mix$model[k])) {
    rows <- which(mix$model[k] == i)
    r <- bets_fit$results[[i]]
    p <- r$particles
    idx <- mix$particle[k[rows]]
    sigma2 <- p$sigma2_scale[idx] / stats::rchisq(length(rows), df = p$sigma2_df)
    # x_T = m_T + \sigma t(R_T) z, z ~ N(0, I)
    z <- matrix(stats::rnorm(dim(p$final_R)[1] * length(rows)), ncol = length(rows))
    states <- t(final_states_rb(p$final_R, p$final_mean, sweep(z, 2, sqrt(sigma2), "*"), idx))
    colnames(states) <- rownames(p$final_mean)
    paths[rows, ] <- ets_future_traj(r$model_components, states, p$theta[idx, , drop = FALSE],
                                     sigma2, h)
  }
  paths
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
