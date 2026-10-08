################################################################################
# INTEGRAND AND MODE (shared by quadrature and AIS)

# Integrand log g(z) = log p(y | \theta(z)) + log p(\theta(z)) + log |d\theta / dz|,
# whose integral is the evidence p(y). Evaluated at the rows of Z in one C++
# call; also returns \theta and ml_res (output of marginal_likelihood_rb(), plus
# the final states as affine functions of eta). With posterior = FALSE, only
# log g: enough for the prior scan and the mode search.
make_log_g_rb <- function(y, model_components, theta_names, ctrl, prior,
                          posterior = TRUE) {
  nu0     <- ctrl$nu0
  phi_min <- ctrl$phi_min
  phi_max <- ctrl$phi_max
  eta0    <- prior$eta0
  V0      <- prior$V0
  psi0    <- prior$psi0
  log_prior_theta_const <- prior$log_prior_theta_const

  L <- length(y)
  m <- stats::frequency(y)
  flags <- model_flags(model_components)
  y_vec <- as.numeric(y)

  function(Z) {
    colnames(Z) <- theta_names
    trans <- transform_unconstrained_to_theta(Z, theta_names, phi_min, phi_max)
    design <- build_design_and_c_batch(
      yR = y_vec,
      trend = flags$trend,
      seas = flags$seas,
      damped = flags$damped,
      m = m,
      paramsR = trans$theta,
      return_final = posterior
    )
    ml_res <- marginal_likelihood_rb(
      XtX_cube = design$XtX,
      Xty_mat  = design$Xty,
      yty_vec  = design$yty,
      eta0     = eta0,
      V0       = V0,
      nu0      = nu0,
      psi0     = psi0,
      L        = L,
      return_posterior = posterior
    )
    if (posterior) {
      ml_res$final_coef  <- design$final_coef
      ml_res$final_const <- design$final_const
    }
    log_g <- as.numeric(ml_res$log_marginal_lik) + log_prior_theta_const + trans$log_jac
    list(log_g = log_g, theta = trans$theta, ml_res = ml_res)
  }
}

# Scan of the prior to start the mode search: a heuristic point, a point with
# small smoothing and n_scan randomised Sobol points, evaluated in one batch.
# Returns the points and their log g, best first. Avoids local modes in the
# unstable region (e.g. monthly AAA).
prior_scan <- function(log_g_fn, theta_names, n_scan = 64) {
  d <- length(theta_names)
  Z <- rbind(heuristic_z_start(theta_names), small_smoothing_z(theta_names))
  if (n_scan > 0) {
    u <- matrix(qrng::sobol(n_scan, d = d, randomize = "digital.shift"), ncol = d)
    Z <- rbind(Z, stats::qlogis(u))
  }
  Z <- Z[rowSums(!is.finite(Z)) == 0, , drop = FALSE]
  log_g <- log_g_fn(Z)$log_g
  log_g[!is.finite(log_g)] <- -Inf
  o <- order(log_g, decreasing = TRUE)
  list(Z = Z[o, , drop = FALSE], log_g = log_g[o])
}

# Heuristic point: centre of each range (z = 0), phi higher (z = 1)
heuristic_z_start <- function(theta_names) {
  as.numeric(theta_names == "phi")
}

# Point with small smoothing parameters
small_smoothing_z <- function(theta_names) {
  z <- c(alpha = stats::qlogis(0.2), beta = stats::qlogis(0.05), phi = 1, gamma = stats::qlogis(0.05))
  unname(z[theta_names])
}

# Mode zhat of log g and Hessian H of -log g there, with Lmat such that
# Lmat t(Lmat) = H^{-1} (eigenvalues of H floored at lambda_floor)
laplace_mode <- function(log_g_fn, z_start,
                         fd_step = 1e-3,
                         lambda_floor = 1e-2,
                         z_bound = 25,
                         penalty = 1e10) {
  d <- length(z_start)
  timing <- list(mode = 0, hessian = 0)

  # -log g, with a finite penalty where log g is not finite
  neg_log_g <- function(Z) {
    f <- -log_g_fn(Z)$log_g
    f[!is.finite(f)] <- penalty
    f
  }

  # Mode (z_bound keeps log_jac finite)
  t0 <- proc.time()[3]
  E_fd <- diag(fd_step, d)
  opt <- stats::optim(
    par = z_start,
    fn  = function(z) neg_log_g(matrix(z, nrow = 1)),
    # Central differences, in one batch
    gr  = function(z) {
      f <- neg_log_g(rbind(sweep(E_fd, 2, z, "+"), sweep(-E_fd, 2, z, "+")))
      (f[seq_len(d)] - f[d + seq_len(d)]) / (2 * fd_step)
    },
    method = "L-BFGS-B",
    lower = -z_bound,
    upper = z_bound
  )
  zhat <- opt$par
  timing$mode <- proc.time()[3] - t0

  # Hessian
  t0 <- proc.time()[3]
  H <- fd_hessian(neg_log_g, zhat, fd_step)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  lambda <- eig$values
  if (any(lambda < lambda_floor)) {
    warning(sprintf(paste0(
      "laplace_mode: Hessian of -log g at the mode is not safely ",
      "positive definite (smallest eigenvalue %.3g); flooring eigenvalues at %g"),
      min(lambda), lambda_floor))
    lambda <- pmax(lambda, lambda_floor)
  }
  Lmat <- eig$vectors %*% diag(lambda^(-1/2), d)
  log_det_H <- sum(log(lambda))  # floored H, as Lmat
  timing$hessian <- proc.time()[3] - t0

  list(zhat = zhat, H = H, Lmat = Lmat, log_det_H = log_det_H, optim = opt, timing = timing)
}

# Hessian of f at z by central differences, all points in one call of f
fd_hessian <- function(f, z, h) {
  d  <- length(z)
  E  <- diag(h, d)
  ij <- which(upper.tri(E), arr.ind = TRUE)   # pairs i < j
  P  <- E[ij[, 1], , drop = FALSE]
  Q  <- E[ij[, 2], , drop = FALSE]
  n_pairs <- nrow(ij)

  # z, z + h e_i, z - h e_i, z + h (+-e_i +-e_j)
  offsets <- rbind(0, E, -E, P + Q, P - Q, -P + Q, -P - Q)
  fv <- f(sweep(offsets, 2, z, "+"))

  f0 <- fv[1]
  fp <- fv[1 + seq_len(d)]
  fm <- fv[1 + d + seq_len(d)]
  H <- diag((fp - 2 * f0 + fm) / h^2, d)
  if (n_pairs > 0) {
    k <- 1 + 2 * d
    f_pp <- fv[k + seq_len(n_pairs)]
    f_pm <- fv[k + n_pairs + seq_len(n_pairs)]
    f_mp <- fv[k + 2 * n_pairs + seq_len(n_pairs)]
    f_mm <- fv[k + 3 * n_pairs + seq_len(n_pairs)]
    H[ij] <- (f_pp - f_pm - f_mp + f_mm) / (4 * h^2)
    H[ij[, 2:1, drop = FALSE]] <- H[ij]
  }
  H
}


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

# Inverse of transform_unconstrained_to_theta(): unconstrained values of the
# smoothing parameters (matrix theta, one column per parameter)
transform_theta_to_unconstrained <- function(theta, phi_min, phi_max) {
  z <- theta
  alpha <- theta[, "alpha"]
  z[, "alpha"] <- stats::qlogis(alpha)
  if ("beta" %in% colnames(theta))  z[, "beta"]  <- stats::qlogis(theta[, "beta"] / alpha)
  if ("gamma" %in% colnames(theta)) z[, "gamma"] <- stats::qlogis(theta[, "gamma"] / (1 - alpha))
  if ("phi" %in% colnames(theta)) {
    z[, "phi"] <- stats::qlogis((theta[, "phi"] - phi_min) / (phi_max - phi_min))
  }
  z
}
