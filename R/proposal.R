##############################################################################
### Transformation Functions (Theta <-> Unconstrained) ###

# Logit and Inverse Logit
logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) 1 / (1 + exp(-x))

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

# Transform Constrained -> Unconstrained
transform_theta_to_unconstrained <- function(theta, phi_min, phi_max) {
  theta_unc <- theta
  eps <- 1e-6 # stability

  if ("alpha" %in% colnames(theta)) {
    val <- pmin(pmax(theta[, "alpha"], eps), 1 - eps)
    theta_unc[, "alpha"] <- logit(val)
  }

  if ("beta" %in% colnames(theta)) {
    # Reconstruct ratio: beta / alpha
    val <- theta[, "beta"] / pmax(theta[, "alpha"], eps)
    val <- pmin(pmax(val, eps), 1 - eps)
    theta_unc[, "beta"] <- logit(val)
  }

  if ("gamma" %in% colnames(theta)) {
    # Reconstruct ratio: gamma / (1 - alpha)
    val <- theta[, "gamma"] / pmax(1 - theta[, "alpha"], eps)
    val <- pmin(pmax(val, eps), 1 - eps)
    theta_unc[, "gamma"] <- logit(val)
  }

  if ("phi" %in% colnames(theta)) {
    val <- (theta[, "phi"] - phi_min) / (phi_max - phi_min)
    val <- pmin(pmax(val, eps), 1 - eps)
    theta_unc[, "phi"] <- logit(val)
  }

  theta_unc
}


##############################################################################
### Joint Proposal Functions ###

# Compute heuristic initialization for parameters of mvt prior for initial states
# Returns list of named vectors (init_mu, init_var)
# names: "l", optionally "b", "s1..sm"
init_eta_params <- function(y, model_components, eta_df,
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
  init_mu <- l
  init_var <- var_l

  if (trend) {
    if (m > 1) {
      b <- (mean(y[(m + 1):(2 * m)]) - mean(y[1:m])) / m
    } else {
      b <- y[2] - y[1]
    }
    var_b <- var_b_mult * stats::var(seas_diffs / n_lags)
    init_mu <- c(init_mu, b)
    init_var <- c(init_var, var_b)
  }

  if (seas && m > 1) {
    s <- rep(0, m - 1)
    y_mat <- matrix(c(y, rep(NA, -L %% m)), nrow = m)
    seas_vars <- apply(y_mat, 2, stats::var, na.rm = TRUE)
    var_s <- rep(var_s_mult * mean(seas_vars, na.rm = TRUE), m - 1)
    init_mu <- c(init_mu, s)
    init_var <- c(init_var, var_s)
  }

  names(init_mu) <- c("l", if (trend) "b", if (seas && m > 1) paste0("s", 1:(m - 1)))
  Sigma <- matrix(0, nrow = length(init_var), ncol = length(init_var))
  diag(Sigma) <- init_var

  list(mus = init_mu, Sigma = Sigma, df = eta_df)
}

# Initialize Joint Proposal
init_joint_params <- function(y, model_components, theta_names, eta_df = 7) {

  # Eta initialization (use heuristic)
  eta_init <- init_eta_params(y, model_components, eta_df = eta_df)
  eta_names <- names(eta_init$mus)

  # Theta initialization (Unconstrained)
  # initialize theta around sensible defaults (e.g. 0 in logit space = 0.5 prob)
  theta_mus <- rep(0, length(theta_names))
  names(theta_mus) <- theta_names

  # If damped trend, phi usually high, set init logit > 0
  if ("phi" %in% theta_names) theta_mus["phi"] <- 1.0

  # Combine
  joint_mus <- c(theta_mus, eta_init$mus)

  # 4. Joint Covariance
  # Block diagonal: Theta part (Identity*scale) + Eta part (Heuristic)
  n_theta <- length(theta_names)
  n_eta <- length(eta_names)

  Sigma_theta <- diag(1, n_theta)

  # Build block diagonal matrix
  Sigma_joint <- matrix(0, nrow = n_theta + n_eta, ncol = n_theta + n_eta)
  Sigma_joint[1:n_theta, 1:n_theta] <- Sigma_theta
  Sigma_joint[(n_theta + 1):(n_theta + n_eta), (n_theta + 1):(n_theta + n_eta)] <- eta_init$Sigma

  rownames(Sigma_joint) <- colnames(Sigma_joint) <- names(joint_mus)

  list(mus = joint_mus, Sigma = Sigma_joint, df = eta_df)
}

# Draw from Joint Proposal
draw_from_joint_proposal <- function(N, proposal_params,
                                     theta_names, phi_min, phi_max,
                                     antithetic = TRUE) {

  df <- proposal_params$df
  mu <- proposal_params$mus
  Sigma <- proposal_params$Sigma
  param_names <- names(mu)

  # 1. Sample Unconstrained Joint Vector
  if (antithetic) {
    base <- mvtnorm::rmvt(n = ceiling(N / 2), sigma = Sigma, df = df, type = "shifted")
    samps_unc <- rbind(base, -base) + matrix(mu, nrow = nrow(base) * 2, ncol = length(mu), byrow = TRUE)
    samps_unc <- samps_unc[1:N, , drop = FALSE]
  } else {
    samps_unc <- mvtnorm::rmvt(n = N, sigma = Sigma, df = df, delta = mu, type = "shifted")
  }
  colnames(samps_unc) <- param_names

  # 2. Compute Proposal Log Density (on unconstrained space)
  log_density_unc <- mvtnorm::dmvt(samps_unc, delta = mu, sigma = Sigma, df = df, log = TRUE)

  # 3. Split into Theta and Eta
  theta_cols <- which(param_names %in% theta_names)
  eta_cols   <- which(!param_names %in% theta_names)

  theta_unc <- samps_unc[, theta_cols, drop = FALSE]
  eta_free  <- samps_unc[, eta_cols,   drop = FALSE]

  # 4. Transform Theta -> Constrained
  trans_res <- transform_unconstrained_to_theta(theta_unc, theta_names, phi_min, phi_max)
  theta_con <- trans_res$theta

  # Adjust log density by Jacobian: log q(theta, eta) = log q(u) - log |J|
  log_density <- log_density_unc - trans_res$log_jac

  # 5. Handle Eta constraints (Seasonal sum to 0)
  # Identify if we have seasonal components s1...sm-1
  s_cols_free <- grep("^s\\d+$", colnames(eta_free), value = TRUE)
  if (length(s_cols_free) > 0) {
    # Calculate last seasonal state
    last_s <- -rowSums(eta_free[, s_cols_free, drop = FALSE])
    eta_con <- cbind(eta_free, last_s)
    # determine m index
    m_idx <- max(as.integer(sub("s", "", s_cols_free))) + 1
    colnames(eta_con)[ncol(eta_con)] <- paste0("s", m_idx)
  } else {
    eta_con <- eta_free
  }

  list(
    theta = theta_con,
    eta = eta_con,
    theta_unc = theta_unc,
    eta_free = eta_free,
    log_density = log_density
  )
}

# Update Joint Proposal
update_joint_proposal <- function(theta_unc, eta_free, w, prev_params,
                                  lr = 0.9, min_var = 1e-6,
                                  lambda_shr = 0.1) {

  # Combine unconstrained samples
  joint_samps <- cbind(theta_unc, eta_free)

  w <- w / sum(w)
  ess <- 1 / sum(w^2)

  # Weighted Mean
  new_mus <- colSums(w * joint_samps)

  # Weighted Covariance
  centered <- sweep(joint_samps, 2, new_mus, "-")
  Sigma_new <- crossprod(centered * sqrt(w))

  # Regularization (Shrinkage + Min Diagonal)
  diag_Sigma <- pmax(diag(Sigma_new), min_var)
  Sigma_new <- (1 - lambda_shr) * Sigma_new
  diag(Sigma_new) <- diag_Sigma

  # Smooth Update
  lamb <- min(lr, ess / (100 + ess))

  out_mus <- lamb * new_mus + (1 - lamb) * prev_params$mus
  out_Sigma <- lamb * Sigma_new + (1 - lamb) * prev_params$Sigma

  list(mus = out_mus, Sigma = out_Sigma, df = prev_params$df)
}

# Vectorized function to generate h-steps ahead trajectories
# for multiple parameter samples with additive ETS models
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

##############################################################################
### NMIG (Normal-Mixture Inverse Gamma) Spike-and-Slab Prior ###

#' Log-density of a spike-and-slab mixture of two zero-mean Normals.
#'
#' @param x     Numeric vector of values.
#' @param v_spike Variance of the spike component (small, e.g. 1e-5).
#' @param v_slab  Variance of the slab component (large, e.g. 10).
#' @param w      Prior mixing weight on the slab (0 < w < 1).
#' @return Numeric vector of log-densities, same length as \code{x}.
log_prior_nmig <- function(x, v_spike, v_slab, w) {
  log_comp_spike <- log(1 - w) + stats::dnorm(x, mean = 0, sd = sqrt(v_spike), log = TRUE)
  log_comp_slab  <- log(w)     + stats::dnorm(x, mean = 0, sd = sqrt(v_slab),  log = TRUE)
  # Numerically stable log-sum-exp over the two components
  max_log <- pmax(log_comp_spike, log_comp_slab)
  max_log + log(exp(log_comp_spike - max_log) + exp(log_comp_slab - max_log))
}
