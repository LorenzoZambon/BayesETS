##############################################################################
### Transformation Functions (Theta <-> Unconstrained) ###

# Logit and Inverse Logit
inv_logit <- stats::plogis   # compiled-C equivalent of 1/(1+exp(-x))

# Transform Unconstrained -> Constrained
# Returns list(theta, log_jac)
# theta_unc: matrix of unconstrained parameters
transform_unconstrained_to_theta <- function(theta_unc, param_names, phi_min, phi_max) {
  N <- nrow(theta_unc)
  theta <- matrix(0, nrow = N, ncol = length(param_names))
  colnames(theta) <- param_names

  # Initialize log-Jacobian
  log_jac <- numeric(N)

  # 1. Alpha: logit(alpha)
  # alpha = inv_logit(x)
  if ("alpha" %in% param_names) {
    p <- inv_logit(theta_unc[, "alpha"])
    theta[, "alpha"] <- p
    log_jac <- log_jac + log(p) + log(1 - p)
  }

  # 2. Beta: logit(beta / alpha)
  if ("beta" %in% param_names) {
    p <- inv_logit(theta_unc[, "beta"])
    # beta = alpha * p
    theta[, "beta"] <- theta[, "alpha"] * p
    # Jacobian adjustment for beta: p * (1-p) * alpha
    log_jac <- log_jac + log(p) + log(1 - p) + log(theta[, "alpha"])
  }

  # 3. Gamma: logit(gamma / (1 - alpha))
  if ("gamma" %in% param_names) {
    p <- inv_logit(theta_unc[, "gamma"])
    # gamma = (1-alpha) * p
    theta[, "gamma"] <- (1 - theta[, "alpha"]) * p
    # Jacobian adjustment: p * (1-p) * (1-alpha)
    log_jac <- log_jac + log(p) + log(1 - p) + log(1 - theta[, "alpha"])
  }

  # 4. Phi: logit((phi - min)/(max - min))
  if ("phi" %in% param_names) {
    p <- inv_logit(theta_unc[, "phi"])
    theta[, "phi"] <- p * (phi_max - phi_min) + phi_min
    log_jac <- log_jac + log(p) + log(1 - p) + log(phi_max - phi_min)
  }

  list(theta = theta, log_jac = log_jac)
}


##############################################################################
### Joint Proposal Functions ###

.VAR_FLOOR <- 1e-6
.FB_MULT <- 1e-4

# Compute heuristic initialization for parameters of mvt prior for initial states
# Returns list of named vectors (init_mu, init_var)
# names: "l", optionally "b", "s1..sm"
init_eta_params <- function(y, model_components, eta_df = NULL,
                            var_l_mult = 1,
                            var_b_mult = 1,
                            var_s_mult = 1) {

  L <- length(y)
  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas  <- (model_components[[3]] == "A")

  n_lags <- ifelse(seas, m, 1)
  seas_diffs <- diff(y, lag = n_lags)
  mse_naive <- mean(seas_diffs^2)

  l <- mean(y[1:m])
  var_l <- var_l_mult * mse_naive
  if (is.na(var_l) || var_l <= 0) var_l <- .VAR_FLOOR  # fallback for constant series

  init_mu <- l
  init_var <- var_l

  if (trend) {
    if (m > 1 && L > m) {
      # Slope between the means of the first two periods; with fewer than 2 m
      # observations, the second window is shifted back to end at L.
      k <- min(m, L - m)
      b <- (mean(y[(k + 1):(k + m)]) - mean(y[1:m])) / k
    } else {
      b <- y[2] - y[1]
    }
    var_b <- var_b_mult * stats::var(seas_diffs / n_lags)
    # ensure b variance is not too small relative to l + fallback if seas_diffs variance is 0
    var_b <- max(var_b, .FB_MULT * var_l, na.rm = T)
    init_mu <- c(init_mu, b)
    init_var <- c(init_var, var_b)
  }

  if (seas && m > 1) {
    s <- rep(0, m - 1)
    y_mat <- matrix(c(y, rep(NA, -L %% m)), nrow = m)
    seas_vars <- apply(y_mat, 2, stats::var, na.rm = TRUE)
    var_s_raw <- var_s_mult * mean(seas_vars, na.rm = TRUE)
    # ensure s variance is not too small relative to l + fallback if seas_diffs variances are 0
    var_s_raw <- max(var_s_raw, .FB_MULT * var_l, na.rm = T)
    var_s <- rep(var_s_raw, m - 1)
    init_mu <- c(init_mu, s)
    init_var <- c(init_var, var_s)
  }

  names(init_mu) <- c("l", if (trend) "b", if (seas && m > 1) paste0("s", 1:(m - 1)))
  Sigma <- matrix(0, nrow = length(init_var), ncol = length(init_var))
  diag(Sigma) <- init_var

  list(mus = init_mu, Sigma = Sigma, df = eta_df)
}

##############################################################################

##############################################################################
### Random-Search Initialization (Sobol QMC) ###

# Evaluate the profile-RSS objective (eta analytically integrated out) for a
# batch of N theta candidates supplied as an N x d unconstrained matrix.
#   failed      : TRUE when the scan could not be completed
sobol_scan_rb <- function(y, model_components, theta_names, phi_min, phi_max,
                          n_sobol, eta0, V0, nu0, psi0, L,
                          log_prior_theta_const, eta_df = 5L) {
  d      <- length(theta_names)
  m      <- stats::frequency(y)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  y_vec  <- as.numeric(y)

  # Generate Sobol points in (0,1)^d with digital-shift randomisation.
  # qlogis maps them to unconstrained space, where each coordinate is
  # Logistic(0,1).  For d >= 2 this proposal is NOT uniform on the admissible
  # ETS region (e.g. gamma = (1 - alpha) * u has density 1 / (1 - alpha)), so
  # the IS weights below use its exact density.  The points are not clipped:
  # clipping piles the tail mass onto the clip boundary, which no density
  # accounts for, and the tails (e.g. alpha near 1) can carry posterior mass.
  # No fixed anchor point is appended either: a deterministic point inside a
  # random sample biases the IS estimate.
  pts_01 <- tryCatch(
    qrng::sobol(n_sobol, d = d, randomize = "digital.shift"),
    error = function(e) NULL
  )
  if (is.null(pts_01)) return(list(prop_params = NULL, ess = 0, failed = TRUE))
  if (d == 1L) pts_01 <- matrix(pts_01, ncol = 1)
  theta_unc_mat <- matrix(stats::qlogis(pts_01), nrow = n_sobol, ncol = d)
  # Drop the (measure-zero) points that land exactly on 0 or 1.
  theta_unc_mat <- theta_unc_mat[rowSums(!is.finite(theta_unc_mat)) == 0, , drop = FALSE]
  colnames(theta_unc_mat) <- theta_names
  n_pts <- nrow(theta_unc_mat)
  if (n_pts == 0L) return(list(prop_params = NULL, ess = 0, failed = TRUE))

  # Transform to constrained space.
  trans         <- transform_unconstrained_to_theta(theta_unc_mat, theta_names,
                                                    phi_min, phi_max)
  theta_con_mat <- trans$theta

  # Build design matrices — single C++ batch call for all n_pts candidates.
  design <- tryCatch(
    build_design_and_c_batch(
      yR = y_vec, trend = trend, seas = seas,
      damped = damped, m = m, paramsR = theta_con_mat
    ),
    error = function(e) NULL
  )
  if (is.null(design)) return(list(prop_params = NULL, ess = 0, failed = TRUE))

  # Marginal likelihood evaluation — same C++ call as a single AIS-RB iteration.
  ml_res <- tryCatch(
    marginal_likelihood_rb(
      XtX_cube = design$XtX, Xty_mat = design$Xty, yty_vec = design$yty,
      eta0 = eta0, V0 = V0, nu0 = nu0, psi0 = psi0, L = L,
      return_posterior = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(ml_res)) return(list(prop_params = NULL, ess = 0, failed = TRUE))

  log_ml <- as.numeric(ml_res$log_marginal_lik)
  log_ml[!is.finite(log_ml)] <- -Inf

  # IS weights w.r.t. Lebesgue measure on theta, as in the AIS loop:
  #   w_i = p(y|θ_i) * p(θ_i) / q(θ_i),   q(θ) = q_z(z) / |dθ/dz|,
  # with q_z the product of Logistic(0,1) densities and log|dθ/dz| = log_jac.
  # mean(w) is then an RQMC estimate of p(y), used by the AIS early exit.
  log_q_unc <- rowSums(stats::dlogis(theta_unc_mat, log = TRUE))
  log_w <- log_ml + log_prior_theta_const + trans$log_jac - log_q_unc
  log_w[!is.finite(log_w)] <- -Inf

  lw_max <- max(log_w[is.finite(log_w)])
  if (!is.finite(lw_max)) {
    warning("sobol_scan_rb: all marginal likelihoods are -Inf; returning failed")
    return(list(prop_params = NULL, ess = 0, failed = TRUE))
  }
  w_raw <- exp(log_w - lw_max)
  w_raw[!is.finite(w_raw)] <- 0
  w   <- w_raw / sum(w_raw)
  ess <- 1 / sum(w^2)
  if (!is.finite(ess)) ess <- 0

  # IS-weighted theta proposal — analogous to one AIS adaptive update step.
  wp <- sobol_weighted_theta_params(theta_unc_mat, w)
  names(wp$mu)                         <- theta_names
  rownames(wp$Sigma) <- colnames(wp$Sigma) <- theta_names
  prop_params <- list(mus = wp$mu, Sigma = wp$Sigma, df = as.integer(eta_df))

  list(
    prop_params = prop_params,
    ess         = ess,
    n_pts       = n_pts,
    theta_unc   = theta_unc_mat,
    theta_con   = theta_con_mat,
    w           = w,
    log_w       = log_w,
    ml_res      = ml_res,
    failed      = FALSE
  )
}

# Random-search single-component joint proposal initialization for AIS.
# "random_search_mean": best single Sobol candidate → mean; heuristic covariance.
# "random_search":      IS-weighted mean *and* covariance from all Sobol candidates,
ets_future_traj <- function(model_components, states, params, sigma2s, h = 10, seed = NULL) {

  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)

  # Model components
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")

  N_samples <- nrow(params)

  # infer m from the provided states
  s_cols <- grep("^s\\d+$", colnames(states), value = TRUE)
  if (length(s_cols) > 1) {
    s_idx <- as.integer(sub("^s", "", s_cols))
    s_cols <- s_cols[order(s_idx)]
  }
  m <- if (seas) length(s_cols) else 1
  if (seas && m == 0) stop("Seasonal model but no seasonal states found.")

  # Extract parameters (already in matrix form)
  alpha <- params[, "alpha"]
  if (trend) {
    beta <- params[, "beta"]
    phi  <- if (damped) params[, "phi"] else 1
  }
  if (seas) gamma <- params[, "gamma"]

  # Initialize states from the RSS_vect output
  l <- states[, "l"]
  if (trend) b <- states[, "b"]
  if (seas)  s <- states[, s_cols, drop = FALSE]

  # Sample forecast errors for all horizons (each sample has its own variance)
  errors <- matrix(stats::rnorm(N_samples * h), nrow = N_samples, ncol = h) * sqrt(sigma2s)

  # Initialize forecast matrix and compute forecasts
  forecasts <- matrix(nrow = N_samples, ncol = h)
  for (i in 1:h) {
    # Compute point forecast for period i
    forecasts[, i] <- l
    if (trend) forecasts[, i] <- forecasts[, i] + phi * b
    if (seas)  forecasts[, i] <- forecasts[, i] + s[, ((i - 1) %% m) + 1]
    # Add forecast error
    forecasts[, i] <- forecasts[, i] + errors[, i]

    # Update states for next period using the realized value (trajectory with errors)
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

log_prior_theta_uniform <- function(theta_samp, phi_min, phi_max) {
  # Initialize with 0 (log density of 1)
  lp <- rep(0, nrow(theta_samp))

  # 1. Beta and Gamma Normalization (Joint Uniform)
  has_beta  <- "beta" %in% colnames(theta_samp)
  has_gamma <- "gamma" %in% colnames(theta_samp)

  if (has_beta && has_gamma) {
    # ETS(A,A,A): Valid volume is 1/6. To integrate to 1, density must be 6.
    lp <- lp + log(6)
  } else if (has_beta || has_gamma) {
    # ETS(A,A,N) or ETS(A,N,A): Valid area is 1/2. Density must be 2.
    lp <- lp + log(2)
  }

  # 2. Phi Normalization (Independent Uniform)
  if ("phi" %in% colnames(theta_samp)) {
    # Density is 1 / (max - min)
    lp <- lp - log(phi_max - phi_min)
  }

  lp
}
