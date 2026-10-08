################################################################################
# PARAMETER TRANSFORMATION

# Smoothing parameters from the unconstrained ones (matrix theta_unc), with the
# log-Jacobian of the transformation: list(theta, log_jac)
transform_unconstrained_to_theta <- function(theta_unc, param_names, phi_min, phi_max) {
  N <- nrow(theta_unc)
  theta <- matrix(0, nrow = N, ncol = length(param_names))
  colnames(theta) <- param_names

  log_jac <- numeric(N)

  # alpha in (0, 1)
  if ("alpha" %in% param_names) {
    p <- stats::plogis(theta_unc[, "alpha"])
    theta[, "alpha"] <- p
    log_jac <- log_jac + log(p) + log(1 - p)
  }

  # beta in (0, alpha)
  if ("beta" %in% param_names) {
    p <- stats::plogis(theta_unc[, "beta"])
    theta[, "beta"] <- theta[, "alpha"] * p
    log_jac <- log_jac + log(p) + log(1 - p) + log(theta[, "alpha"])
  }

  # gamma in (0, 1 - alpha)
  if ("gamma" %in% param_names) {
    p <- stats::plogis(theta_unc[, "gamma"])
    theta[, "gamma"] <- (1 - theta[, "alpha"]) * p
    log_jac <- log_jac + log(p) + log(1 - p) + log(1 - theta[, "alpha"])
  }

  # phi in (phi_min, phi_max)
  if ("phi" %in% param_names) {
    p <- stats::plogis(theta_unc[, "phi"])
    theta[, "phi"] <- p * (phi_max - phi_min) + phi_min
    log_jac <- log_jac + log(p) + log(1 - p) + log(phi_max - phi_min)
  }

  list(theta = theta, log_jac = log_jac)
}


################################################################################
# PRIOR OF THE INITIAL STATES

.VAR_FLOOR <- 1e-6   # floor of the prior variance of l
.FB_MULT <- 1e-4     # floor of the prior variances of b and s, relative to l

# Heuristic prior mean and (diagonal) covariance of the initial states
# (l, b, s1, ..., s_{m-1}); var_*_mult scale the variances (for prior tuning)
init_eta_params <- function(y, model_components,
                            var_l_mult = 1,
                            var_b_mult = 1,
                            var_s_mult = 1) {

  L <- length(y)
  m <- stats::frequency(y)
  flags <- model_flags(model_components)
  trend <- flags$trend
  seas  <- flags$seas

  n_lags <- ifelse(seas, m, 1)
  seas_diffs <- diff(y, lag = n_lags)
  mse_naive <- mean(seas_diffs^2)

  l <- mean(y[1:m])
  var_l <- var_l_mult * mse_naive
  if (is.na(var_l) || var_l <= 0) var_l <- .VAR_FLOOR

  init_mu <- l
  init_var <- var_l

  if (trend) {
    if (m > 1 && L > m) {
      # Slope between the means of the first two periods (the second one
      # shifted back if L < 2m)
      k <- min(m, L - m)
      b <- (mean(y[(k + 1):(k + m)]) - mean(y[1:m])) / k
    } else {
      b <- y[2] - y[1]
    }
    var_b <- var_b_mult * stats::var(seas_diffs / n_lags)
    var_b <- max(var_b, .FB_MULT * var_l, na.rm = T)
    init_mu <- c(init_mu, b)
    init_var <- c(init_var, var_b)
  }

  if (seas && m > 1) {
    # Prior variance: mean variance within a period
    s <- rep(0, m - 1)
    y_mat <- matrix(c(y, rep(NA, -L %% m)), nrow = m)
    seas_vars <- apply(y_mat, 2, stats::var, na.rm = TRUE)
    var_s_raw <- var_s_mult * mean(seas_vars, na.rm = TRUE)
    var_s_raw <- max(var_s_raw, .FB_MULT * var_l, na.rm = T)
    var_s <- rep(var_s_raw, m - 1)
    init_mu <- c(init_mu, s)
    init_var <- c(init_var, var_s)
  }

  names(init_mu) <- c("l", if (trend) "b", if (seas && m > 1) paste0("s", 1:(m - 1)))
  Sigma <- matrix(0, nrow = length(init_var), ncol = length(init_var))
  diag(Sigma) <- init_var

  list(mus = init_mu, Sigma = Sigma)
}

################################################################################
# FORECASTS AND PRIOR

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

# Log-density of the uniform prior of the smoothing parameters (a constant)
log_prior_theta_uniform <- function(theta_names, phi_min, phi_max) {
  # (alpha, beta, gamma): uniform on a region of volume 1/6 (AAA) or 1/2
  # (AAN, ANA)
  has_beta  <- "beta" %in% theta_names
  has_gamma <- "gamma" %in% theta_names
  lp <- if (has_beta && has_gamma) log(6) else if (has_beta || has_gamma) log(2) else 0

  # phi: uniform on (phi_min, phi_max)
  if ("phi" %in% theta_names) lp <- lp - log(phi_max - phi_min)
  lp
}
