#' Fitted values and residuals of a BETS model
#'
#' One-step-ahead fitted values and residuals on the training series.
#'
#' @param object A fitted object from [bets()].
#' @param type `"mean"` (default): posterior mean of the fitted values, over the
#'   models (with their weights) and over the posterior of their parameters
#'   (computed exactly, without sampling); the residuals are `y` minus this
#'   mean. `"draws"`: one row per posterior draw, from the posterior of the
#'   model combination.
#' @param n_draws Number of rows for `type = "draws"` (default 1000).
#' @param ... Unused.
#'
#' @return With `type = "mean"`, a `ts` with the time index of `y`. With
#'   `type = "draws"`, a matrix with `n_draws` rows and one column per
#'   observation.
#'
#' @examples
#' set.seed(1)
#' fit <- bets(USAccDeaths)
#'
#' # One-step-ahead fitted values (posterior mean) and residuals
#' fitted(fit)
#' residuals(fit)
#'
#' # Draws of the fitted values (from the posterior) and residuals
#' fitted(fit, type = "draws")[1:5,1:12]   # only show first 5 draws and 12 observations
#' residuals(fit, type = "draws")[1:5,1:12]
#'
#' @export
fitted.bets <- function(object, type = c("mean", "draws"), n_draws = 1000, ...) {
  type <- match.arg(type)
  y <- as.numeric(object$y)
  res <- one_step_residuals(object, type, n_draws)
  if (type == "mean") {
    stats::ts(y - res, start = stats::tsp(object$y)[1], frequency = stats::frequency(object$y))
  } else {
    sweep(-res, 2, y, "+")
  }
}

#' @rdname fitted.bets
#' @export
residuals.bets <- function(object, type = c("mean", "draws"), n_draws = 1000, ...) {
  type <- match.arg(type)
  res <- one_step_residuals(object, type, n_draws)
  if (type == "mean") {
    stats::ts(res, start = stats::tsp(object$y)[1], frequency = stats::frequency(object$y))
  } else {
    res
  }
}

# One-step residuals, from the particles of the models. "mean": given theta, the
# residuals are affine in the initial states eta, so their posterior mean is the
# residuals at the posterior mean mu_n of eta; averaged over particles and models.
# "draws": n_draws rows, each from a particle drawn by its weight (model weight x
# particle weight), with \sigma^2 and eta from their posteriors.
one_step_residuals <- function(object, type, n_draws) {
  if (type == "draws" && !is_count(n_draws)) {
    stop("n_draws must be a single positive integer")
  }
  results <- object$fit$results
  w <- object$fit$model_weights
  L <- length(object$y)
  # Models with particles (not failed) and positive weight
  models <- which(w > 0 & !vapply(results, function(r) is.null(r$particles), logical(1)))

  if (type == "draws") {
    n_k <- vapply(results[models], function(r) length(r$particles$w), integer(1))
    model <- rep(models, n_k)
    particle <- sequence(n_k)
    pw <- unlist(lapply(models, function(i) w[i] * results[[i]]$particles$w))
    k <- sample.int(length(pw), as.integer(n_draws), replace = TRUE, prob = pw)
    out <- matrix(0, n_draws, L)
  } else {
    out <- numeric(L)
  }

  for (i in models) {
    r <- results[[i]]
    # Constant series: zero residuals
    if (r$integration == "constant series") next
    p <- r$particles
    pp <- particle_posteriors(object, i)
    ml <- pp$ml_res
    if (type == "mean") {
      idx <- seq_along(p$w)
      eta <- ml$mu_n
    } else {
      rows <- which(model[k] == i)
      if (length(rows) == 0) next
      idx <- particle[k[rows]]
      sigma2 <- p$sigma2_scale[idx] / stats::rchisq(length(rows), df = p$sigma2_df)
      # eta = mu_n + \sigma t(Rn) z, z ~ N(0, I)
      z <- matrix(stats::rnorm(nrow(ml$mu_n) * length(rows)), ncol = length(rows))
      eta <- final_states_rb(ml$Rn, ml$mu_n, sweep(z, 2, sqrt(sigma2), "*"), idx)
    }
    E <- residuals_at(pp$y, r$model_components, object$period, eta, p$theta[idx, , drop = FALSE])
    if (type == "mean") out <- out + w[i] * colSums(p$w * E) else out[rows, ] <- E
  }
  out
}

# One-step residuals (one row per column of eta) of an additive ETS model on the
# series y, from the initial states eta (C++ order) with smoothing parameters theta
residuals_at <- function(y, model_components, m, eta, theta) {
  flags <- model_flags(model_components)
  ets_residuals(as.numeric(y), flags$trend, flags$seas, flags$damped, if (flags$seas) m else 1L,
                eta_to_r_order(eta, model_components, m), theta, return_residuals = TRUE)$residuals
}
