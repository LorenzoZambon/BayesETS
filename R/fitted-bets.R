#' Fitted values and residuals of a BETS model
#'
#' One-step-ahead fitted values and residuals on the training series.
#'
#' @param object A fitted object from [bets()].
#' @param type `"mean"` (default): posterior mean of the fitted values, averaged
#'   over the posterior draws of each model and over the models with their
#'   weights; the residuals are `y` minus this mean. `"draws"`: one row per
#'   posterior draw, with the draws of each model in proportion to its weight.
#' @param n_draws Number of rows for `type = "draws"` (default 1000).
#' @param ... Unused.
#'
#' @return With `type = "mean"`, a `ts` with the time index of `y`. With
#'   `type = "draws"`, a matrix with `n_draws` rows and one column per
#'   observation.
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

# One-step residuals of the posterior draws: their weighted mean over draws and
# models ("mean"), or n_draws rows of draws, resampled by model weight ("draws")
one_step_residuals <- function(object, type, n_draws) {
  if (type == "draws" && (!is.numeric(n_draws) || length(n_draws) != 1 || n_draws < 1)) {
    stop("n_draws must be a single positive number")
  }
  y <- as.numeric(object$y)
  results <- object$fit$results
  w <- object$fit$model_weights
  n_k <- if (type == "draws") draws_per_model(w, as.integer(n_draws))

  out <- lapply(seq_along(results), function(k) {
    r <- results[[k]]
    if (w[k] == 0 || (type == "draws" && n_k[k] == 0)) return(NULL)
    idx <- if (type == "draws") {
      sample(nrow(r$thetas), n_k[k], replace = TRUE)
    } else {
      seq_len(nrow(r$thetas))
    }
    flags <- model_flags(r$model_components)
    m <- max(1L, sum(grepl("^s\\d+$", colnames(r$etas))))   # seasonal period
    E <- RSS_vect_arma(y, flags$trend, flags$seas, flags$damped, m,
                       r$etas[idx, , drop = FALSE], r$thetas[idx, , drop = FALSE],
                       return_residuals = TRUE)$residuals
    if (type == "mean") w[k] * colMeans(E) else E
  })
  out <- Filter(Negate(is.null), out)

  if (type == "mean") Reduce(`+`, out) else do.call(rbind, out)
}
