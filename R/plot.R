#' Plot the posterior of a BETS model
#'
#' One figure per model: the posterior densities of the smoothing parameters and
#' of `sigma`, one panel each, with the credible interval shaded and the
#' posterior mean as a dashed line. With `states = TRUE`, a second figure per
#' model shows the initial level and slope, and the seasonal states as violins,
#' ordered and labelled by season (month or quarter) when possible.
#'
#' The densities of the smoothing parameters are weighted kernel density
#' estimates from the posterior particles, reflected at the bounds of the
#' parameters. Those of `sigma` and of the initial states are exact mixtures
#' over the particles. The initial states are normalised as in [summary.bets()].
#'
#' @param x A fitted object from [bets()].
#' @param models Labels of the models to plot (e.g. `"AAA"`). Default: the
#'   models with weight of at least 0.001.
#' @param level Level of the credible intervals, in percent (default 95).
#' @param states Logical: also plot the initial states (default `FALSE`).
#' @param ask Logical: ask before each new figure, if there is more than one
#'   (default: in interactive sessions).
#' @param ... Unused.
#'
#' @return Invisibly, a data frame with the plotted densities: `model`,
#'   `parameter`, `x`, `density` and `in_interval`.
#'
#' @examples
#' set.seed(1)
#' fit <- bets(USAccDeaths)
#' plot(fit, models = "ANA")
#' plot(fit, models = "ANA", states = TRUE)
#'
#' @export
plot.bets <- function(x, models = NULL, level = 95, states = FALSE,
                      ask = grDevices::dev.interactive(), ...) {
  if (!(is.numeric(level) && length(level) == 1 && is.finite(level) && level > 0 && level < 100))
    stop("level must be a single number between 0 and 100")
  if (!(is.logical(states) && length(states) == 1 && !is.na(states)))
    stop("states must be TRUE or FALSE")

  labels <- vapply(x$model_components, ets_label, character(1))
  results <- x$fit$results
  w <- x$fit$model_weights
  if (is.null(models)) {
    idx <- shown_models(w)
  } else {
    if (!is.character(models) || anyNA(match(models, labels)))
      stop("models must be labels of the fitted models: ", paste(labels, collapse = ", "))
    idx <- match(models, labels)
  }
  failed <- vapply(results[idx], function(r) is.null(r$particles), logical(1))
  if (any(failed))
    stop("No posterior for the failed model(s): ", paste(labels[idx][failed], collapse = ", "))

  probs <- c(100 - level, 100 + level) / 200
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (isTRUE(ask) && length(idx) * (1 + states) > 1) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }

  curves <- list()
  for (i in idx) {
    title <- sprintf("%s: weight %.3f", labels[i], w[i])
    p <- results[[i]]$particles
    s <- posterior_summary(p, probs)

    # Smoothing parameters and sigma, in 2 columns
    graphics::par(mfrow = c(ceiling(nrow(s) / 2), 2), oma = c(0, 0, 2, 0), mar = c(3, 3, 2.5, 1))
    for (par_name in rownames(s)) {
      d <- if (par_name == "sigma") {
        sigma_density(p)
      } else {
        bounds <- if (par_name == "phi") c(x$control$phi_min, x$control$phi_max) else c(0, 1)
        parameter_density(p$theta[, par_name], p$w, bounds)
      }
      interval <- unlist(s[par_name, 3:4])
      density_panel(d, s[par_name, "mean"], interval, parse(text = par_name))
      curves[[length(curves) + 1]] <- curve_data(labels[i], par_name, d, interval)
    }
    graphics::mtext(title, outer = TRUE, font = 2, line = 0.5)

    if (states) {
      curves <- c(curves, plot_states(x, i, probs, paste(title, "- initial states")))
    }
  }
  invisible(do.call(rbind, curves))
}

# Figure of the initial states of model i: level and slope as densities, the
# seasonal states as violins ordered by season. Returns the curves (data frames).
plot_states <- function(object, i, probs, title) {
  p <- object$fit$results[[i]]$particles
  st <- initial_states_particles(object, i)
  states <- colnames(st$loc)
  dens <- lapply(states, function(j) t_mixture_density(st$loc[, j], st$scale[, j], p$sigma2_df, p$w))
  summ <- t(vapply(states, function(j) {
    t_mixture_summary(st$loc[, j], st$scale[, j], p$sigma2_df, p$w, probs)
  }, numeric(4)))
  names(dens) <- states

  seasonal <- grep("^s", states, value = TRUE)
  other <- setdiff(states, seasonal)
  top <- c(seq_along(other), if (length(other) == 1) 0)
  if (length(seasonal) > 0) {
    graphics::layout(rbind(top, length(other) + 1), heights = c(1, 1.3))
  } else {
    graphics::layout(rbind(top))
  }
  graphics::par(oma = c(0, 0, 2, 0), mar = c(3, 3, 2.5, 1))
  names_other <- c(l = "initial level", b = "initial slope")
  for (j in other) density_panel(dens[[j]], summ[j, 1], summ[j, 3:4], names_other[[j]])
  if (length(seasonal) > 0) {
    seasons <- season_order(object$y, length(seasonal))
    violin_panel(dens[seasons$states], summ[seasons$states, 1], summ[seasons$states, 3:4, drop = FALSE],
                 seasons$labels, "initial seasonal states")
  }
  graphics::mtext(title, outer = TRUE, font = 2, line = 0.5)

  lapply(states, function(j) {
    curve_data(ets_label(object$fit$results[[i]]$model_components), j, dens[[j]], summ[j, 3:4])
  })
}

# Seasonal states s1, ..., s_m ordered by season, with labels: s_m applies to
# the first observation, s_{m-1} to the second, ...; month or quarter names if
# the time index of y allows it, otherwise s1, ..., s_m
season_order <- function(y, m) {
  states <- paste0("s", m - seq_len(m) + 1)   # states of observations 1, ..., m
  if (stats::frequency(y) != m) return(list(states = paste0("s", seq_len(m)), labels = paste0("s", seq_len(m))))
  pos <- stats::cycle(y)[seq_len(m)]
  o <- order(pos)
  labels <- if (m == 12) month.abb[pos] else if (m == 4) paste0("Q", pos) else paste0("s", m - seq_len(m) + 1)
  list(states = states[o], labels = labels[o])
}

# Weighted kernel density of a parameter in the interval `bounds`, from its
# particles, reflected at the bounds (so that no mass leaks outside them).
# Bandwidth by Silverman's rule, with the effective number of particles.
parameter_density <- function(x, w, bounds, n = 512) {
  keep <- w > 0
  x <- x[keep]
  w <- w[keep] / sum(w[keep])
  sd_x <- sqrt(sum(w * (x - sum(w * x))^2))
  iqr_x <- diff(weighted_quantile(x, w, c(0.25, 0.75))) / 1.34
  spread <- if (iqr_x > 0) min(sd_x, iqr_x) else sd_x
  bw <- 0.9 * spread * sum(w^2)^(1 / 5)                   # n_eff = 1 / sum(w^2)
  k <- stats::density(c(x, 2 * bounds[1] - x, 2 * bounds[2] - x), weights = rep(w, 3) / 3,
                      bw = bw, n = n, from = bounds[1], to = bounds[2])
  trim_density(k$x, 3 * k$y)
}

# Exact density of sigma: mixture over the particles of the density of
# sqrt(scale / \chi^2_df)
sigma_density <- function(p, n = 512) {
  q <- sigma_summary(p$sigma2_scale, p$sigma2_df, p$w, c(0.0005, 0.9995))[3:4]
  s <- seq(q[1], q[2], length.out = n)
  d <- vapply(s, function(v) {
    sum(p$w * stats::dchisq(p$sigma2_scale / v^2, p$sigma2_df) * 2 * p$sigma2_scale / v^3)
  }, numeric(1))
  list(x = s, density = d)
}

# Exact density of a mixture of location-scale t (NULL for a point mass)
t_mixture_density <- function(loc, scale, df, w, n = 512) {
  if (all(scale == 0)) return(NULL)
  q <- t_mixture_summary(loc, scale, df, w, c(0.0005, 0.9995))[3:4]
  x <- seq(q[1], q[2], length.out = n)
  d <- vapply(x, function(v) sum(w * stats::dt((v - loc) / scale, df) / scale), numeric(1))
  list(x = x, density = d)
}

# The part of a density curve between its 0.05% and 99.95% quantiles
trim_density <- function(x, d) {
  cdf <- c(0, cumsum(diff(x) * (d[-1] + d[-length(d)]) / 2))
  cdf <- cdf / cdf[length(cdf)]
  keep <- seq(max(1, sum(cdf < 0.0005)), min(length(x), sum(cdf <= 0.9995) + 1))
  list(x = x[keep], density = d[keep])
}

# Density panel: curve, shaded credible interval, dashed line at the mean (a
# vertical line for a point mass)
density_panel <- function(d, mean, interval, main, col = "steelblue") {
  if (is.null(d)) {
    graphics::plot(mean, 1, type = "h", lwd = 2, col = col, main = main, xlab = "", ylab = "", yaxt = "n")
    return(invisible())
  }
  graphics::plot(d$x, d$density, type = "n", main = main, xlab = "", ylab = "",
                 ylim = c(0, 1.05 * max(d$density)), yaxs = "i",
                 cex.main = if (is.expression(main)) 1.6 else 1.2)
  xs <- c(interval[1], d$x[d$x > interval[1] & d$x < interval[2]], interval[2])
  ys <- stats::approx(d$x, d$density, xs, rule = 2)$y
  graphics::polygon(c(xs[1], xs, xs[length(xs)]), c(0, ys, 0),
                    col = grDevices::adjustcolor(col, 0.35), border = NA)
  graphics::lines(d$x, d$density, col = col, lwd = 2)
  graphics::abline(v = mean, lty = 2)
}

# Violins of several densities (list of curves, drawn vertically), with the
# credible intervals shaded and the means as points
violin_panel <- function(dens, means, intervals, labels, main, col = "steelblue") {
  ylim <- range(unlist(lapply(dens, `[[`, "x")))
  dmax <- max(unlist(lapply(dens, `[[`, "density")))
  graphics::plot(NA, xlim = c(0.5, length(dens) + 0.5), ylim = ylim, xaxt = "n",
                 xlab = "", ylab = "", main = main)
  graphics::axis(1, at = seq_along(dens), labels = labels)
  graphics::abline(h = 0, col = "grey60")
  for (j in seq_along(dens)) {
    d <- dens[[j]]
    hw <- 0.45 * d$density / dmax
    graphics::polygon(c(j - hw, rev(j + hw)), c(d$x, rev(d$x)),
                      col = grDevices::adjustcolor(col, 0.12), border = NA)
    inside <- d$x >= intervals[j, 1] & d$x <= intervals[j, 2]
    graphics::polygon(c(j - hw[inside], rev(j + hw[inside])), c(d$x[inside], rev(d$x[inside])),
                      col = grDevices::adjustcolor(col, 0.35), border = NA)
    graphics::polygon(c(j - hw, rev(j + hw)), c(d$x, rev(d$x)), col = NA, border = col)
    graphics::points(j, means[j], pch = 19, cex = 0.7)
  }
}

# Data frame of a plotted curve
curve_data <- function(model, parameter, d, interval) {
  if (is.null(d)) {
    return(data.frame(model = model, parameter = parameter, x = interval[[1]], density = NA_real_,
                      in_interval = TRUE))
  }
  data.frame(model = model, parameter = parameter, x = d$x, density = d$density,
             in_interval = d$x >= interval[[1]] & d$x <= interval[[2]])
}
