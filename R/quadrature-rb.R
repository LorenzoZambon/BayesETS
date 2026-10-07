################################################################################
# ADAPTIVE GAUSS-HERMITE QUADRATURE
#
# Alternative to AIS: the smoothing parameters are integrated on a
# Gauss-Hermite grid at the posterior mode, scaled by the inverse Hessian (in
# the unconstrained space). The initial states and \sigma^2 are integrated
# analytically, as in AIS.

# Quadrature over the smoothing parameters, from z_start: log evidence and
# weighted nodes (theta, w, ml_res); failure: message, or NULL
integrate_quadrature <- function(log_g_fn, z_start, theta_names, ctrl,
                                 log_g_mode = log_g_fn) {
  n_quad <- resolve_by_d(ctrl$n_quad, length(theta_names))
  quad <- adaptive_gh_quadrature(log_g_fn, z_start, n_quad, log_g_mode = log_g_mode)

  if (isTRUE(ctrl$verbose >= 2)) {
    cat(sprintf("\n\nRB Gauss-Hermite quadrature: %d nodes (%d per dimension)\n",
                nrow(quad$G), n_quad))
    cat(sprintf("\n log evidence = %.4f\n", quad$log_evidence))
  }

  list(
    log_evidence = quad$log_evidence,
    theta = quad$node_eval$theta,
    w = quad$w,
    ml_res = quad$node_eval$ml_res,
    ess = NA_real_,
    n_iter = NA_integer_,
    proposal = NULL,
    timing = quad$timing,
    failure = if (!is.finite(quad$log_evidence)) {
      "Quadrature: log g is not finite at any node; the model gets zero weight"
    }
  )
}


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


################################################################################
# QUADRATURE HELPERS

# Adaptive Gauss-Hermite quadrature of g = exp(log g) on R^d, with n_quad nodes
# per dimension. log_g_fn(Z) returns a list with the values log_g at the rows
# of Z (other elements are kept in node_eval); z_start starts the mode search,
# done with log_g_mode (same values, possibly cheaper). Nodes z_i = zhat + Lmat x_i, and
#   \int g(z) dz = |Lmat| E[g(zhat + Lmat X) / \phi_d(X)],  X ~ N(0, I_d)
adaptive_gh_quadrature <- function(log_g_fn, z_start, n_quad,
                                   fd_step = 1e-3,
                                   lambda_floor = 1e-2,
                                   z_bound = 25,
                                   penalty = 1e10,
                                   log_g_mode = log_g_fn) {
  d <- length(z_start)

  # Mode and Hessian
  lap <- laplace_mode(log_g_mode, z_start, fd_step = fd_step, lambda_floor = lambda_floor,
                      z_bound = z_bound, penalty = penalty)
  zhat <- lap$zhat
  Lmat <- lap$Lmat
  log_det_H <- lap$log_det_H
  timing <- c(lap$timing, list(nodes = 0))

  # Grid
  gh  <- gauss_hermite_prob(n_quad)
  idx <- as.matrix(expand.grid(rep(list(seq_len(n_quad)), d)))
  G   <- matrix(gh$x[idx], ncol = d)
  log_wq <- rowSums(matrix(log(gh$w[idx]), ncol = d))  # log-weights of the nodes
  Z <- sweep(G %*% t(Lmat), 2, zhat, "+")

  # log g at all nodes, in one batch
  t0 <- proc.time()[3]
  node_eval <- log_g_fn(Z)
  log_g <- node_eval$log_g
  log_g[!is.finite(log_g)] <- -Inf
  timing$nodes <- proc.time()[3] - t0

  # log of |Lmat| g(Z_i) / \phi_d(G_i)
  corr <- log_g + 0.5 * rowSums(G^2) - 0.5 * log_det_H + (d / 2) * log(2 * pi)
  log_terms <- log_wq + corr
  lt_max <- max(log_terms)
  if (is.finite(lt_max)) {
    w <- exp(log_terms - lt_max)
    log_evidence <- lt_max + log(sum(w))
    w <- w / sum(w)
  } else {
    w <- NULL
    log_evidence <- -Inf
  }

  list(
    log_evidence = log_evidence,
    w = w,
    Z = Z,
    G = G,
    node_eval = node_eval,
    zhat = zhat,
    H = lap$H,
    log_det_H = log_det_H,
    optim = lap$optim,
    timing = timing
  )
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

# Gauss-Hermite rule for E[f(X)], X ~ N(0, 1) (statmod gives the rule for the
# weight function exp(-x^2))
gauss_hermite_prob <- function(k) {
  gq <- statmod::gauss.quad(k, kind = "hermite")
  x <- sqrt(2) * gq$nodes
  w <- gq$weights / sqrt(pi)
  stopifnot("Gauss-Hermite weights must sum to 1" = abs(sum(w) - 1) < 1e-10)
  list(x = x, w = w)
}
