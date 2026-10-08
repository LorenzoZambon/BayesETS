################################################################################
# PRIOR (shared by quadrature and AIS)

# Prior of the initial states and \sigma^2, and log-density of the uniform
# prior of the smoothing parameters. eta0 and V0 are in the C++ order.
init_rb_prior <- function(y, model_components, theta_names, ctrl) {
  psi0          <- ctrl$psi0
  phi_min       <- ctrl$phi_min
  phi_max       <- ctrl$phi_max
  c_inflate_eta <- ctrl$c_inflate_eta

  m <- stats::frequency(y)
  flags <- model_flags(model_components)
  trend <- flags$trend
  seas <- flags$seas

  eta_init <- init_eta_params(y, model_components)  # l, [b,] [s1, ..., s_{m-1}]

  # eta also includes the last seasonal state: (l0, [b0,] [s1, ..., s_m])
  n_eta <- n_states(model_components, m)

  # Prior mean of s_m: -(s1 + ... + s_{m-1})
  eta0_free <- eta_init$mus
  if (seas && m > 1) {
    s_free_names <- grep("^s\\d+$", names(eta0_free), value = TRUE)
    last_s <- -sum(eta0_free[s_free_names])
    eta0_r_order <- c(eta0_free, last_s)
    names(eta0_r_order)[length(eta0_r_order)] <- paste0("s", m)
  } else {
    eta0_r_order <- eta0_free
  }

  # eta | \sigma^2 ~ N(eta0, \sigma^2 V0), with V0 = Sigma / E[\sigma^2]: prior
  # covariance Sigma at the prior mean of \sigma^2
  Sigma_heuristic_free <- eta_init$Sigma * c_inflate_eta
  prior_mean_sigma2 <- psi0 / (ctrl$nu0 - 2)

  # Independent priors on all m seasonal states, without the sum-to-zero
  # constraint: the prior fixes the shift (l0 + c, s - c) that the data cannot
  # identify, and keeps V0 invertible
  if (seas && m > 1) {
    n_free <- length(eta0_free)
    n_full <- n_eta
    s_indices_free <- grep("^s\\d+$", names(eta0_free))
    Sigma_full <- matrix(0, n_full, n_full)
    Sigma_full[1:n_free, 1:n_free] <- Sigma_heuristic_free
    # s_m: same prior variance as the other seasonal states
    Sigma_full[n_full, n_full] <- mean(diag(Sigma_heuristic_free)[s_indices_free])
  } else {
    Sigma_full <- Sigma_heuristic_free
  }

  # C++ order of the seasonal states: (s_m, ..., s1)
  if (seas && m > 1) {
    s_offset_r <- (1 + (if (trend) 1L else 0L))  # number of non-seasonal states
    non_s_idx <- seq_len(s_offset_r)
    s_idx_r <- (s_offset_r + 1):(s_offset_r + m)  # s1, ..., s_m
    s_idx_cpp <- rev(s_idx_r)                     # s_m, ..., s1
    reorder <- c(non_s_idx, s_idx_cpp)

    eta0_cpp <- as.numeric(eta0_r_order[reorder])
    V0 <- Sigma_full[reorder, reorder, drop = FALSE] / prior_mean_sigma2
  } else {
    eta0_cpp <- as.numeric(eta0_r_order)
    V0 <- Sigma_full / prior_mean_sigma2
  }

  log_prior_theta_const <- log_prior_theta_uniform(theta_names, phi_min, phi_max)

  list(
    eta0 = eta0_cpp,
    V0 = V0,
    psi0 = psi0,
    log_prior_theta_const = log_prior_theta_const
  )
}

# Default psi0: variance of the naive or (if smaller) of the seasonal naive residuals,
# floored to stay positive for deterministic series; heuristic for constant series.
default_psi0 <- function(y) {
  if (all(y == y[1])) return((0.2 * (if (y[1] != 0) abs(y[1]) else 1))^2)
  m <- stats::frequency(y)
  psi0 <- stats::var(diff(y))
  if (m > 1) psi0 <- min(psi0, stats::var(diff(y, lag = m)), na.rm = TRUE)
  max(psi0, 1e-8 * stats::var(y))
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
