##############################################################################
### Rao-Blackwellized Adaptive Gauss-Hermite Quadrature ###
###
### Alternative to adaptive_is_rb(): theta is integrated out deterministically
### on a tensor Gauss-Hermite grid centred at the mode of the integrand and
### scaled by its inverse Hessian, in the unconstrained coordinates z of
### transform_unconstrained_to_theta().  eta and sigma^2 are integrated out
### analytically exactly as in adaptive_is_rb().
##############################################################################

quadrature_rb <- function(y, model_components, ctrl,
                          return_pointwise = FALSE) {
  N_final_raw <- ctrl$N_final
  n_quad_raw  <- ctrl$n_quad    # scalar or length-4 vector; resolved per d below
  nu0         <- ctrl$nu0
  phi_min     <- ctrl$phi_min
  phi_max     <- ctrl$phi_max
  verbose     <- ctrl$verbose

  L <- length(y)
  m <- stats::frequency(y)
  trend <- (model_components[[2]] == "A")
  seas <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha", if (trend) c("beta", if (damped) "phi"), if (seas) "gamma")
  n_theta <- length(theta_names)

  # Resolve dimension-dependent scalars now that d = n_theta is known.
  N_final <- resolve_by_d(N_final_raw, n_theta)
  n_quad  <- resolve_by_d(n_quad_raw,  n_theta)

  # ---- Set up eta prior ----
  prior <- init_rb_prior(y, model_components, theta_names, ctrl)
  eta0  <- prior$eta0
  V0    <- prior$V0
  psi0  <- prior$psi0
  log_prior_theta_const <- prior$log_prior_theta_const

  y_vec <- as.numeric(y)

  # log g(z) = log p(y | theta(z)) + log p(theta(z)) + log |d theta / d z|,
  # for every row of Z in one batched C++ call.  The integral of g over z is
  # the model evidence p(y).
  log_g_rb <- function(Z) {
    colnames(Z) <- theta_names
    trans <- transform_unconstrained_to_theta(Z, theta_names, phi_min, phi_max)
    design <- build_design_and_c_batch(
      yR = y_vec,
      trend = trend,
      seas = seas,
      damped = damped,
      m = m,
      paramsR = trans$theta
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
      return_posterior = TRUE
    )
    log_g <- as.numeric(ml_res$log_marginal_lik) + log_prior_theta_const + trans$log_jac
    list(log_g = log_g, theta = trans$theta, ml_res = ml_res)
  }

  # ---- Start of the mode search ----
  # Heuristic start, as in init_joint_params(): z = 0 (centre of each range),
  # except phi at z = 1 (damping usually high).  No Sobol scan: on M3 it gave
  # the same log evidence (to 1e-4) at a large share of the cost.
  z_start <- as.numeric(theta_names == "phi")

  # ---- Quadrature ----
  quad <- adaptive_gh_quadrature(log_g_rb, z_start, n_quad)
  timing <- c(quad$timing, list(post = 0))

  if (verbose >= 2) {
    cat(sprintf("\n\nRB Gauss-Hermite quadrature: %d nodes (%d per dimension)\n",
                nrow(quad$G), n_quad))
    cat(sprintf("\n log evidence = %.4f\n", quad$log_evidence))
  }

  # No finite node: this model will receive zero weight in BMA/stacking.
  if (!is.finite(quad$log_evidence)) {
    warning("quadrature_rb: log g is not finite at any quadrature node")
    return(list(
      thetas = NULL,
      etas = NULL,
      states = NULL,
      sigma2s = NULL,
      ess = NA_real_,
      n_iter = NA_integer_,
      prop_params = NULL,
      log_evidence = -Inf,
      log_lik_pointwise = if (return_pointwise) matrix(-1e300, nrow = N_final, ncol = L) else NULL,  # use -1e300 instead of -Inf to avoid NaN in logsumexp
      timing = timing
    ))
  }

  # ---- Posterior Reconstruction ----
  t0 <- proc.time()[3]
  post <- draw_rb_posterior(y, model_components,
                            theta_particles = quad$node_eval$theta,
                            w = quad$w,
                            ml_res = quad$node_eval$ml_res,
                            N_final = N_final, nu0 = nu0,
                            return_pointwise = return_pointwise)
  timing$post <- proc.time()[3] - t0

  list(
    thetas = post$thetas,
    etas = post$etas,
    states = post$states,
    sigma2s = post$sigma2s,
    ess = NA_real_,
    n_iter = NA_integer_,
    prop_params = NULL,
    log_evidence = quad$log_evidence,
    log_lik_pointwise = post$log_lik_pointwise,
    timing = timing
  )
}


##############################################################################
### Quadrature Helpers ###

# Adaptive Gauss-Hermite quadrature of g = exp(log g) over R^d.
#   log_g_fn : function(Z) of an n x d matrix of points (one per row) returning
#              a list whose element log_g holds the n values of log g.  Any other
#              elements are passed through for the nodes (see node_eval).
#   z_start  : starting point of the mode search.
#   n_quad   : number of one-dimensional nodes; the grid has n_quad^d nodes.
# The grid is z_i = zhat + Lmat %*% x_i with Lmat %*% t(Lmat) = H^{-1}, where
# zhat is the mode of log g and H the Hessian of -log g there.  Then
#   integral g(z) dz = |Lmat| * E[ g(zhat + Lmat X) / phi_d(X) ],  X ~ N(0, I_d),
# which is exact for any Lmat; the expectation is taken by Gauss-Hermite.
adaptive_gh_quadrature <- function(log_g_fn, z_start, n_quad,
                                   fd_step = 1e-3,
                                   lambda_floor = 1e-2,
                                   z_bound = 25,
                                   penalty = 1e10) {
  d <- length(z_start)
  timing <- list(mode = 0, hessian = 0, nodes = 0)

  # Finite-valued -log g for the optimiser and the difference stencils: optim
  # may wander into regions where the marginal likelihood is -Inf or NaN.
  neg_log_g <- function(Z) {
    f <- -log_g_fn(Z)$log_g
    f[!is.finite(f)] <- penalty
    f
  }

  # ---- Step 1: Mode ----
  # z_bound keeps log_jac finite (inv_logit saturates to 0/1 beyond |z| ~ 37).
  t0 <- proc.time()[3]
  E_fd <- diag(fd_step, d)
  opt <- stats::optim(
    par = z_start,
    fn  = function(z) neg_log_g(matrix(z, nrow = 1)),
    # Central differences, all 2d points in one batch.
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

  # ---- Step 2: Curvature ----
  t0 <- proc.time()[3]
  H <- fd_hessian(neg_log_g, zhat, fd_step)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  lambda <- eig$values
  if (any(lambda < lambda_floor)) {
    warning(sprintf(paste0(
      "adaptive_gh_quadrature: Hessian of -log g at the mode is not safely ",
      "positive definite (smallest eigenvalue %.3g); flooring eigenvalues at %g"),
      min(lambda), lambda_floor))
    lambda <- pmax(lambda, lambda_floor)
  }
  Lmat <- eig$vectors %*% diag(lambda^(-1/2), d)
  log_det_H <- sum(log(lambda))  # of the floored H, consistent with Lmat
  timing$hessian <- proc.time()[3] - t0

  # ---- Steps 3-4: Grid ----
  gh  <- gauss_hermite_prob(n_quad)
  idx <- as.matrix(expand.grid(rep(list(seq_len(n_quad)), d)))
  G   <- matrix(gh$x[idx], ncol = d)
  log_wq <- rowSums(matrix(log(gh$w[idx]), ncol = d))  # log of row products of w
  Z <- sweep(G %*% t(Lmat), 2, zhat, "+")

  # ---- Step 5: Evaluate log g at all nodes in one batch ----
  t0 <- proc.time()[3]
  node_eval <- log_g_fn(Z)
  log_g <- node_eval$log_g
  log_g[!is.finite(log_g)] <- -Inf
  timing$nodes <- proc.time()[3] - t0

  # ---- Step 6: Combine ----
  # corr_i = log[ |Lmat| * g(Z_i) / phi_d(G_i) ]
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
    H = H,
    log_det_H = log_det_H,
    optim = opt,
    timing = timing
  )
}

# Hessian of f at z by central differences.  f takes an n x d matrix of points
# and returns n values; the whole stencil (2 d^2 + 1 points) is one call.
fd_hessian <- function(f, z, h) {
  d  <- length(z)
  E  <- diag(h, d)
  ij <- which(upper.tri(E), arr.ind = TRUE)   # pairs i < j (none if d = 1)
  P  <- E[ij[, 1], , drop = FALSE]
  Q  <- E[ij[, 2], , drop = FALSE]
  n_pairs <- nrow(ij)

  # Rows: z, z + h e_i, z - h e_i, then z + h(+-e_i +-e_j) for each pair.
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

# One-dimensional Gauss-Hermite rule in the probabilists' convention:
# sum(w * f(x)) ~ E[f(X)], X ~ N(0, 1).  statmod::gauss.quad() returns the
# physicists' rule (weight function exp(-x^2), weights summing to sqrt(pi)).
gauss_hermite_prob <- function(k) {
  gq <- statmod::gauss.quad(k, kind = "hermite")
  x <- sqrt(2) * gq$nodes
  w <- gq$weights / sqrt(pi)
  stopifnot("Gauss-Hermite weights must sum to 1" = abs(sum(w) - 1) < 1e-10)
  list(x = x, w = w)
}
