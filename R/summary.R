#' Posterior summaries of a BETS model
#'
#' Posterior mean, standard deviation and equal-tailed credible interval of the
#' smoothing parameters and of `sigma` (standard deviation of the errors), and
#' optionally of the initial states, for each model. They are computed from the
#' weighted posterior particles of each model (quadrature nodes or AIS draws).
#'
#' The initial states are the level `l`, the slope `b` and the seasonal states
#' `s1`, ..., `sm` (`sm` is the one of the first observation), 
#' normalised to sum to zero.
#'
#' @param object A fitted object from [bets()].
#' @param level Level of the credible intervals, in percent (default 95).
#' @param states Logical: also summarise the initial states (default `FALSE`).
#' @param ... Unused.
#'
#' @return An object of class `summary.bets`: a list with the `call`, the
#'   `level`, a data frame `models` (model, weight, log evidence), a list
#'   `parameters` with one data frame per model (rows: parameters and sigma;
#'   columns: mean, sd and interval bounds; `NULL` for failed models), with
#'   `states = TRUE` a list `states` of the same form, and, with BMA, the
#'   posterior `probabilities` of trend, damped trend and seasonality.
#'   Printing it shows the models with weight of at least 0.001.
#'
#' @examples
#' set.seed(1)
#' fit <- bets(USAccDeaths)
#' summary(fit)
#'
#' # 80% intervals of one model, and its initial states
#' s <- summary(fit, level = 80, states = TRUE)
#' s$parameters$ANA
#' s$states$ANA
#'
#' @export
summary.bets <- function(object, level = 95, states = FALSE, ...) {
  if (!(is.numeric(level) && length(level) == 1 && is.finite(level) && level > 0 && level < 100))
    stop("level must be a single number between 0 and 100")
  if (!(is.logical(states) && length(states) == 1 && !is.na(states)))
    stop("states must be TRUE or FALSE")

  labels <- vapply(object$model_components, ets_label, character(1))
  results <- object$fit$results
  w <- object$fit$model_weights
  # f(i) for the models with particles (NULL for failed models)
  by_model <- function(f) {
    stats::setNames(lapply(seq_along(results), function(i) {
      if (!is.null(results[[i]]$particles)) f(i)
    }), labels)
  }

  structure(
    list(
      call = object$call,
      level = level,
      models = data.frame(model = labels, weight = w,
                          log_evidence = vapply(results, `[[`, numeric(1), "log_evidence")),
      parameters = by_model(function(i) posterior_summary(results[[i]]$particles, level / 100)),
      states = if (states) by_model(function(i) states_summary(object, i, level / 100)),
      probabilities = if (object$fit$combination == "bma") {
        component_probabilities(object$model_components, w)
      }
    ),
    class = "summary.bets"
  )
}

#' @export
print.summary.bets <- function(x, ...) {
  cat("Posterior summaries of a BETS model fit\n")
  cat(sprintf("  call: %s\n", paste(deparse(x$call, width.cutoff = 500L), collapse = " ")))

  m <- x$models
  shown <- shown_models(m$weight)
  for (i in shown) {
    le <- m$log_evidence[i]
    cat(sprintf("\n  %s: weight %.3f%s\n", m$model[i], m$weight[i],
                if (is.finite(le)) sprintf(", log evidence %.2f", le) else ""))
    # Smoothing parameters with 3 decimals; sigma and the states (scale of y)
    # with 4 significant digits
    s <- x$parameters[[i]]
    par <- rownames(s) != "sigma"
    tab <- matrix("", nrow(s), ncol(s), dimnames = dimnames(s))
    tab[par, ] <- sprintf("%.3f", as.matrix(s[par, ]))
    tab[!par, ] <- format_signif(as.matrix(s[!par, ]))
    cat_table(cbind(" " = rownames(s), tab))
    if (!is.null(x$states)) {
      cat("  initial states:\n")
      st <- as.matrix(x$states[[i]])
      cat_table(cbind(" " = rownames(st), matrix(format_signif(st), nrow(st), dimnames = dimnames(st))))
    }
  }
  cat_negligible(m$model, shown)
  cat_probabilities(x$probabilities)
  invisible(x)
}

# Posterior mean, sd and equal-tailed interval (level in (0, 1)) of the
# smoothing parameters and of sigma, from the weighted particles of a model
posterior_summary <- function(particles, level = 0.95) {
  theta <- particles$theta
  w <- particles$w / sum(particles$w)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  rows <- lapply(colnames(theta), function(p) {
    m <- sum(w * theta[, p])
    c(m, sqrt(sum(w * (theta[, p] - m)^2)), weighted_quantile(theta[, p], w, probs))
  })
  rows[[length(rows) + 1]] <- sigma_summary(particles$sigma2_scale, particles$sigma2_df, w, probs)
  out <- as.data.frame(do.call(rbind, rows))
  dimnames(out) <- list(c(colnames(theta), "sigma"), summary_colnames(probs))
  out
}

# Posterior mean, sd and equal-tailed interval of the normalised initial states
# of model i: mixture over the particles of their t posteriors
states_summary <- function(object, i, level = 0.95) {
  p <- object$fit$results[[i]]$particles
  w <- p$w / sum(p$w)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  st <- initial_states_particles(object, i)
  out <- t(vapply(colnames(st$loc), function(j) {
    t_mixture_summary(st$loc[, j], st$scale[, j], p$sigma2_df, w, probs)
  }, numeric(4)))
  out <- as.data.frame(out)
  colnames(out) <- summary_colnames(probs)
  out
}

# Location and scale (particles x states) of the t posteriors (df = sigma2_df)
# of the normalised initial states of model i, at its particles. Not stored in
# the fit, to keep bets() fast: the posterior is evaluated again at the
# particles, on the series as fitted (period and centring of fit_bets_models())
initial_states_particles <- function(object, i) {
  r <- object$fit$results[[i]]
  p <- r$particles
  y <- stats::ts(as.numeric(object$y), frequency = object$period)
  if (r$integration == "constant series") {
    const <- function(v) matrix(v, nrow = length(p$w), ncol = 1, dimnames = list(NULL, "l"))
    return(list(loc = const(y[1]), scale = const(0)))
  }

  ctrl <- object$control
  ctrl$psi0 <- object$psi0
  shift <- centring_shift(y)
  y <- y - shift
  mc <- r$model_components
  tn <- colnames(p$theta)
  log_g <- make_log_g_rb(y, mc, tn, ctrl, init_rb_prior(y, mc, tn, ctrl))
  ml <- log_g(transform_theta_to_unconstrained(p$theta, ctrl$phi_min, ctrl$phi_max))$ml_res

  st <- initial_states_t(ml$mu_n, ml$Rn, as.numeric(ml$posterior_scale) / p$sigma2_df,
                         initial_states_map(mc, object$period))
  st$loc[, "l"] <- st$loc[, "l"] + shift
  st
}

# Map from the initial states in the C++ order (l, [b,] s_m, ..., s1) to the R
# order (l, [b,] s1, ..., s_m), normalised as in forecast::ets(): the seasonal
# states sum to zero and the level absorbs their mean (the data identify only
# l + s_j and the differences of the s_j)
initial_states_map <- function(model_components, m) {
  flags <- model_flags(model_components)
  cpp_names <- c("l", if (flags$trend) "b", if (flags$seas) paste0("s", rev(seq_len(m))))
  r_names   <- c("l", if (flags$trend) "b", if (flags$seas) paste0("s", seq_len(m)))
  k <- length(r_names)
  A <- diag(k)
  if (flags$seas) {
    s <- grep("^s", r_names)
    A[s, s] <- A[s, s] - 1 / m
    A[1, s] <- 1 / m
  }
  A <- A %*% diag(k)[match(r_names, cpp_names), , drop = FALSE]
  dimnames(A) <- list(r_names, cpp_names)
  A
}

# Marginal t posteriors of A eta at each particle (column of mu_n, slice of Rn):
# location A mu_n and scale sqrt(s2 diag(A V_n A')), with V_n = t(Rn) Rn and
# s2 = posterior_scale / df. Returns list(loc, scale), particles x states.
initial_states_t <- function(mu_n, Rn, s2, A) {
  n <- nrow(mu_n)
  N <- ncol(mu_n)
  # Rn A' of all particles in one product (rows of RA: row of Rn, particle)
  RA <- matrix(aperm(Rn, c(1, 3, 2)), ncol = n) %*% t(A)
  v <- colSums(array(RA^2, dim = c(n, N, nrow(A))))
  loc <- t(A %*% mu_n)
  colnames(loc) <- colnames(v) <- rownames(A)
  list(loc = loc, scale = sqrt(v * s2))
}

summary_colnames <- function(probs) c("mean", "sd", paste0(signif(100 * probs, 4), "%"))

# Mean, sd and quantiles of a mixture (weights w) of location-scale t
# distributions with df degrees of freedom (point masses if all scales are 0)
t_mixture_summary <- function(loc, scale, df, w, probs) {
  m <- sum(w * loc)
  v <- sum(w * (scale^2 * df / (df - 2) + loc^2)) - m^2
  q <- if (all(scale == 0)) weighted_quantile(loc, w, probs) else vapply(probs, function(p) {
    # The quantile of the mixture lies between those of the components
    range_p <- range(loc + scale * stats::qt(p, df))
    if (diff(range_p) <= 1e-10 * max(abs(range_p))) return(range_p[1])
    cdf <- function(x) sum(w * stats::pt((x - loc) / scale, df))
    stats::uniroot(function(x) cdf(x) - p, range_p, tol = 1e-10 * max(abs(range_p)))$root
  }, numeric(1))
  c(m, sqrt(max(v, 0)), q)
}

# 4 significant digits, keeping trailing zeros
format_signif <- function(x) sub("\\.$", "", formatC(x, digits = 4, format = "fg", flag = "#"))

# Quantiles of a discrete distribution (values x, weights w), by linear
# interpolation of the CDF at the midpoints of the weights
weighted_quantile <- function(x, w, probs) {
  keep <- w > 0
  o <- order(x[keep])
  x <- x[keep][o]
  w <- w[keep][o] / sum(w[keep])
  stats::approx(cumsum(w) - w / 2, x, xout = probs, rule = 2, ties = "ordered")$y
}

# Mean, sd and quantiles of sigma, whose posterior is a mixture over the
# particles (weights w) of \sigma^2 = scale / \chi^2_df
sigma_summary <- function(scale, df, w, probs) {
  # E[\sigma | theta] = sqrt(scale) E[1 / \chi_df], E[\sigma^2 | theta] = scale / (df - 2)
  m  <- sum(w * sqrt(scale)) * exp(lgamma((df - 1) / 2) - lgamma(df / 2)) / sqrt(2)
  m2 <- sum(w * scale) / (df - 2)
  cdf <- function(s) sum(w * stats::pchisq(scale / s^2, df, lower.tail = FALSE))
  q <- vapply(probs, function(p) {
    # The quantile of the mixture lies between those of the components
    range_p <- sqrt(range(scale) / stats::qchisq(1 - p, df))
    if (diff(range_p) <= 1e-10 * range_p[2]) return(range_p[1])
    stats::uniroot(function(s) cdf(s) - p, range_p, tol = 1e-10 * range_p[2])$root
  }, numeric(1))
  c(m, sqrt(max(m2 - m^2, 0)), q)
}

# BMA: posterior probability of each component (trend, damped trend,
# seasonality) present in some of the models but not all; NULL if none
component_probabilities <- function(model_components, weights) {
  has <- list(trend          = vapply(model_components, `[[`, character(1), 2) != "N",
              `damped trend` = vapply(model_components, `[[`, character(1), 4) == "TRUE",
              seasonality    = vapply(model_components, `[[`, character(1), 3) != "N")
  has <- Filter(function(h) any(h) && !all(h), has)
  if (length(has) > 0) vapply(has, function(h) sum(weights[h]), numeric(1))
}

# Models shown by print() and summary(): weight >= 0.001, by decreasing weight
shown_models <- function(weights) {
  o <- order(weights, decreasing = TRUE)
  o[weights[o] >= 0.001]
}

# The models not shown, in one line
cat_negligible <- function(labels, shown) {
  hidden <- setdiff(seq_along(labels), shown)
  if (length(hidden) > 0) {
    cat(sprintf("  (weight < 0.001: %s)\n", paste(labels[hidden], collapse = ", ")))
  }
}

# "Posterior probability of trend: 0.75, ..." (nothing if probs is NULL)
cat_probabilities <- function(probs) {
  if (length(probs) > 0) {
    cat(sprintf("\n  Posterior probability of %s\n",
                paste(sprintf("%s: %.2f", names(probs), probs), collapse = ", ")))
  }
}
