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
    cat(sprintf("RB Gauss-Hermite quadrature: %d nodes (%d per dimension)",
                nrow(quad$G), n_quad))
    cat(sprintf("\nlog evidence = %.4f\n", quad$log_evidence))
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

# Gauss-Hermite rule for E[f(X)], X ~ N(0, 1) (statmod gives the rule for the
# weight function exp(-x^2))
gauss_hermite_prob <- function(k) {
  gq <- statmod::gauss.quad(k, kind = "hermite")
  x <- sqrt(2) * gq$nodes
  w <- gq$weights / sqrt(pi)
  stopifnot("Gauss-Hermite weights must sum to 1" = abs(sum(w) - 1) < 1e-10)
  list(x = x, w = w)
}
