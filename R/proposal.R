##############################################################################
### Transformation Functions (Theta <-> Unconstrained) ###

# Logit and Inverse Logit
logit <- function(p) log(p / (1 - p))
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

# Compute profile RSS from normal-equation sufficient statistics.
# Some seasonal / short-series designs are singular or very ill-conditioned,
# so an exact solve(XtX, Xty) can fail even when the candidate is usable.
# For initialization we only need a stable relative score, so we fall back to a
# lightly regularized solve before declaring the candidate degenerate.
profile_rss_from_suff_stats <- function(XtX, Xty, yty) {
  XtX <- as.matrix(XtX)
  Xty <- as.numeric(Xty)
  n_eta <- nrow(XtX)
  if (length(Xty) != n_eta || !is.finite(yty)) return(1e15)

  eta_hat <- tryCatch(solve(XtX, Xty), error = function(e) NULL)

  if (is.null(eta_hat)) {
    diag_scale <- mean(diag(XtX))
    if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1.0

    ridge_grid <- diag_scale * c(1e-10, 1e-8, 1e-6)
    for (ridge in ridge_grid) {
      eta_hat <- tryCatch(
        solve(XtX + diag(ridge, n_eta), Xty),
        error = function(e) NULL
      )
      if (!is.null(eta_hat)) break
    }
  }

  if (is.null(eta_hat)) return(1e15)

  rss <- yty - 2 * sum(Xty * eta_hat) + drop(crossprod(eta_hat, XtX %*% eta_hat))
  if (!is.finite(rss)) return(1e15)

  max(rss, 0)
}

##############################################################################
### MLE-Based Initialization ###

# Objective function for Nelder-Mead: given theta in unconstrained space,
# analytically solve for the optimal eta (OLS) and return the resulting RSS.
# Uses build_design_and_c_batch to construct design matrices, avoiding explicit
# eta sampling and enabling fast evaluation.
mle_rss_objective <- function(theta_unc, theta_names, phi_min, phi_max,
                              y_vec, trend, seas, damped, m) {
  theta_mat <- matrix(theta_unc, nrow = 1)
  colnames(theta_mat) <- theta_names
  trans     <- transform_unconstrained_to_theta(theta_mat, theta_names, phi_min, phi_max)
  theta_con <- trans$theta

  design <- tryCatch(
    build_design_and_c_batch(
      yR = y_vec, trend = trend, seas = seas,
      damped = damped, m = m, paramsR = theta_con
    ),
    error = function(e) NULL
  )
  if (is.null(design)) return(1e15)

  XtX <- design$XtX[, , 1]
  Xty <- design$Xty[, 1]
  yty <- design$yty[1]

  profile_rss_from_suff_stats(XtX, Xty, yty)
}

# Find the approximate MLE for theta by running Nelder-Mead, analytically
# integrating out eta at each function evaluation.
# Returns the unconstrained theta at the optimum (named vector).
find_theta_mle <- function(y, model_components, phi_min, phi_max,
                           mle_tol, mle_maxit) {
  m      <- stats::frequency(y)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha",
                   if (trend) c("beta", if (damped) "phi"),
                   if (seas) "gamma")
  y_vec <- as.numeric(y)

  theta_unc_start <- rep(0, length(theta_names))
  names(theta_unc_start) <- theta_names
  if ("phi" %in% theta_names) theta_unc_start["phi"] <- 1.0

  obj <- function(x) {
    mle_rss_objective(x, theta_names, phi_min, phi_max, y_vec, trend, seas, damped, m)
  }

  # Clamp unconstrained result to (-unc_max, unc_max) so that constrained
  # parameters stay away from 0/1 boundaries.  Beyond this range the logit
  # Jacobian approaches -Inf, making the IS log-density +Inf.
  # inv_logit(4) ≈ 0.982, inv_logit(-4) ≈ 0.018 — well within the valid space.
  unc_max <- 10

  if (length(theta_unc_start) == 1L) {
    # optimize() is reliable and exact for 1D; Nelder-Mead is not.
    result <- tryCatch(
      stats::optimize(f = obj, interval = c(-unc_max, unc_max), tol = mle_tol),
      error = function(e) NULL
    )
    if (is.null(result)) {
      warning("MLE initialization via optimize() failed; falling back to heuristic")
      return(theta_unc_start)
    }
    theta_unc_hat <- result$minimum
    names(theta_unc_hat) <- theta_names
    return(theta_unc_hat)
  }

  result <- tryCatch(
    stats::optim(
      par     = theta_unc_start,
      fn      = obj,
      method  = "Nelder-Mead",
      control = list(reltol = mle_tol, maxit = mle_maxit)
    ),
    error = function(e) NULL
  )

  if (is.null(result)) {
    warning("MLE initialization via Nelder-Mead failed; falling back to heuristic")
    return(theta_unc_start)
  }

  theta_unc_hat <- pmin(pmax(result$par, -unc_max), unc_max)
  names(theta_unc_hat) <- theta_names
  theta_unc_hat
}

# MLE-based joint proposal initialization.
# Centers the theta part of the proposal at the approximate MLE; keeps the
# heuristic initialization for eta and the covariance structure unchanged.
init_joint_params_mle <- function(y, model_components, theta_names, eta_df,
                                  phi_min, phi_max, mle_tol, mle_maxit) {
  base_params        <- init_joint_params(y, model_components, theta_names, eta_df)
  theta_unc_hat      <- find_theta_mle(y, model_components, phi_min, phi_max,
                                       mle_tol, mle_maxit)
  base_params$mus[theta_names] <- theta_unc_hat
  base_params
}

##############################################################################
### Random-Search Initialization (Sobol QMC) ###

# Evaluate the profile-RSS objective (eta analytically integrated out) for a
# batch of N theta candidates supplied as an N x d unconstrained matrix.
# All N design matrices are built in a single C++ call; the OLS system is
# then solved for each candidate in a small R loop (at most 512 tiny solves).
# Returns a numeric vector of length N.
eval_rss_batch <- function(theta_unc_mat, theta_names, phi_min, phi_max,
                           y_vec, trend, seas, damped, m) {
  n <- nrow(theta_unc_mat)

  trans     <- transform_unconstrained_to_theta(theta_unc_mat, theta_names,
                                                phi_min, phi_max)
  theta_con <- trans$theta

  # Helper: solve the OLS system for one design-matrix slot and return RSS.
  # rss is clamped to 0 from below: a tiny negative value is floating-point
  # cancellation indicating a near-perfect fit, not a degenerate case.
  solve_rss_one <- function(XtX_i, Xty_i, yty_i) {
    profile_rss_from_suff_stats(XtX_i, Xty_i, yty_i)
  }

  # Fast path: evaluate all N candidates in a single C++ batch call.
  design <- tryCatch(
    build_design_and_c_batch(
      yR = y_vec, trend = trend, seas = seas,
      damped = damped, m = m, paramsR = theta_con
    ),
    error = function(e) NULL
  )

  if (!is.null(design)) {
    return(vapply(seq_len(n), function(i) {
      solve_rss_one(design$XtX[, , i], design$Xty[, i], design$yty[i])
    }, numeric(1)))
  }

  # Fallback: the batch call threw (one bad candidate can cause C++ to raise).
  # Evaluate each candidate individually so that isolated failures only produce
  # a single 1e15 entry rather than discarding the entire batch.
  vapply(seq_len(n), function(i) {
    d_i <- tryCatch(
      build_design_and_c_batch(
        yR = y_vec, trend = trend, seas = seas,
        damped = damped, m = m,
        paramsR = theta_con[i, , drop = FALSE]
      ),
      error = function(e) NULL
    )
    if (is.null(d_i)) return(1e15)
    solve_rss_one(d_i$XtX[, , 1], d_i$Xty[, 1], d_i$yty[1])
  }, numeric(1))
}

# Evaluate the Rao-Blackwellized marginal log-likelihood for a batch of theta
# candidates. This uses the same design-matrix statistics as the RB sampler,
# so the random-search objective matches the downstream target more closely
# than plain profile RSS.
eval_rb_log_ml_batch <- function(theta_unc_mat, theta_names, phi_min, phi_max,
                                 y_vec, trend, seas, damped, m, rb_scoring) {
  n <- nrow(theta_unc_mat)

  trans     <- transform_unconstrained_to_theta(theta_unc_mat, theta_names,
                                                phi_min, phi_max)
  theta_con <- trans$theta

  eval_one <- function(XtX_i, Xty_i, yty_i) {
    ml_res <- tryCatch(
      marginal_likelihood_rb(
        XtX_cube = array(XtX_i, dim = c(nrow(XtX_i), ncol(XtX_i), 1L)),
        Xty_mat  = matrix(Xty_i, ncol = 1L),
        yty_vec  = yty_i,
        eta0     = rb_scoring$eta0,
        V0       = rb_scoring$V0,
        nu0      = rb_scoring$nu0,
        psi0     = rb_scoring$psi0,
        L        = rb_scoring$L
      ),
      error = function(e) NULL
    )
    if (is.null(ml_res)) return(-Inf)
    log_ml <- as.numeric(ml_res$log_marginal_lik)[1]
    if (!is.finite(log_ml)) return(-Inf)
    log_ml
  }

  design <- tryCatch(
    build_design_and_c_batch(
      yR = y_vec, trend = trend, seas = seas,
      damped = damped, m = m, paramsR = theta_con
    ),
    error = function(e) NULL
  )

  if (!is.null(design)) {
    ml_res <- tryCatch(
      marginal_likelihood_rb(
        XtX_cube = design$XtX,
        Xty_mat  = design$Xty,
        yty_vec  = design$yty,
        eta0     = rb_scoring$eta0,
        V0       = rb_scoring$V0,
        nu0      = rb_scoring$nu0,
        psi0     = rb_scoring$psi0,
        L        = rb_scoring$L
      ),
      error = function(e) NULL
    )

    if (!is.null(ml_res)) {
      log_ml <- as.numeric(ml_res$log_marginal_lik)
      log_ml[!is.finite(log_ml)] <- -Inf
      return(log_ml)
    }
  }

  vapply(seq_len(n), function(i) {
    d_i <- tryCatch(
      build_design_and_c_batch(
        yR = y_vec, trend = trend, seas = seas,
        damped = damped, m = m,
        paramsR = theta_con[i, , drop = FALSE]
      ),
      error = function(e) NULL
    )
    if (is.null(d_i)) return(-Inf)
    eval_one(d_i$XtX[, , 1], d_i$Xty[, 1], d_i$yty[1])
  }, numeric(1))
}

default_n_sobol <- function(d) {
  as.integer(2L ^ (d + 3L))   # 16, 32, 64, 128 for d = 1..4
}

# Find the best theta starting point(s) by evaluating a Sobol low-discrepancy
# sequence over the constrained parameter space.  Sobol points are generated
# in [0,1]^d and transformed via qlogis() so that they map uniformly over the
# admissible constrained region:
#   alpha          ~ U(0,1)
#   beta/alpha     ~ U(0,1)
#   gamma/(1-alpha)~ U(0,1)
#   phi            ~ U(phi_min, phi_max)
# This is superior to uniform sampling in unconstrained space, which would
# produce a U-shaped (boundary-heavy) distribution after inv_logit.
#
# n_sobol: NULL → auto-select as 2^(d+5) = 64/128/256/512 for d=1..4.
#          Powers of 2 are optimal for Sobol sequences.
# top_k  : number of best candidates to return (rows of returned matrix).
#          top_k = 1 returns the single best point as a named vector.
#
# Note: qrng is used instead of randtoolbox because it is purpose-built for
# quasi-random generation, has a simpler/stateless API, and the digital-shift
# randomization gives proper RQMC variance estimates.
find_theta_random_search <- function(y, model_components, phi_min, phi_max,
                                     n_sobol = NULL, top_k = 1L,
                                     score = c("rss", "rb_marglik"),
                                     rb_scoring = NULL) {
  score <- match.arg(score)
  m      <- stats::frequency(y)
  trend  <- (model_components[[2]] == "A")
  seas   <- (model_components[[3]] == "A")
  damped <- (model_components[[4]] == "TRUE")
  theta_names <- c("alpha",
                   if (trend) c("beta", if (damped) "phi"),
                   if (seas) "gamma")
  y_vec <- as.numeric(y)
  d     <- length(theta_names)

  if (is.null(n_sobol)) n_sobol <- default_n_sobol(d)

  if (score == "rb_marglik" && is.null(rb_scoring)) {
    stop("rb_scoring must be provided when score = 'rb_marglik'")
  }

  heuristic_row <- {
    hr <- rep(0, d)
    names(hr) <- theta_names
    if ("phi" %in% theta_names) hr["phi"] <- 1.0
    hr
  }

  fallback_mat <- matrix(rep(heuristic_row, top_k), nrow = top_k, byrow = TRUE)
  colnames(fallback_mat) <- theta_names

  # Sobol points in (0,1)^d with digital-shift randomization.
  # qlogis maps them uniformly over the admissible constrained region.
  # The heuristic center is appended as an extra candidate so random_search
  # never performs worse than simply evaluating the heuristic point alone.
  pts_01 <- tryCatch(
    qrng::sobol(n_sobol, d = d, randomize = "digital.shift"),
    error = function(e) NULL
  )
  if (is.null(pts_01)) {
    warning("random_search: qrng::sobol() failed; falling back to heuristic center")
    return(if (top_k == 1L) heuristic_row else fallback_mat)
  }

  # Map Sobol points to unconstrained space via qlogis, then clamp to
  # [-unc_max, unc_max].  The clamp is essential: without it, Sobol points
  # near 0 or 1 produce qlogis values of ±Inf (or large magnitudes like ±10),
  # which cause C++ overflow/NaN inside build_design_and_c_batch and bring down
  # the entire batch.  With unc_max = 4, constrained coverage is
  # inv_logit([-4,4]) ≈ [0.018, 0.982] — uniform over the interior of the
  # admissible ETS parameter region.
  unc_max <- 4
  if (d == 1L) pts_01 <- matrix(pts_01, ncol = 1)
  theta_unc_mat <- matrix(
    pmin(pmax(stats::qlogis(pts_01), -unc_max), unc_max),
    nrow = n_sobol, ncol = d
  )
  colnames(theta_unc_mat) <- theta_names
  theta_unc_mat <- rbind(theta_unc_mat, heuristic_row)
  rownames(theta_unc_mat) <- NULL

  n_candidates <- nrow(theta_unc_mat)
  top_k <- min(top_k, n_candidates)

  score_vec <- switch(score,
    rss = eval_rss_batch(theta_unc_mat, theta_names, phi_min, phi_max,
                         y_vec, trend, seas, damped, m),
    rb_marglik = eval_rb_log_ml_batch(theta_unc_mat, theta_names, phi_min, phi_max,
                                      y_vec, trend, seas, damped, m,
                                      rb_scoring = rb_scoring)
  )

  all_bad <- switch(score,
    rss = all(score_vec >= 1e15),
    rb_marglik = all(!is.finite(score_vec) | score_vec == -Inf)
  )
  if (all_bad) {
    warning("random_search: all candidates invalid; falling back to heuristic center")
    return(if (top_k == 1L) heuristic_row else fallback_mat)
  }

  best_idx <- switch(score,
    rss = order(score_vec)[seq_len(top_k)],
    rb_marglik = order(score_vec, decreasing = TRUE)[seq_len(top_k)]
  )
  result   <- theta_unc_mat[best_idx, , drop = FALSE]
  rownames(result) <- NULL

  if (top_k == 1L) result[1L, ] else result
}

# Random-search single-component joint proposal initialization (for AIS).
# Centers the theta part of the proposal at the best Sobol candidate; keeps
# the heuristic covariance structure unchanged.
init_joint_params_random_search <- function(y, model_components, theta_names,
                                            eta_df, phi_min, phi_max, n_sobol,
                                            score = "rss", rb_scoring = NULL) {
  base_params   <- init_joint_params(y, model_components, theta_names, eta_df)
  theta_unc_hat <- find_theta_random_search(y, model_components,
                                            phi_min, phi_max,
                                            n_sobol = n_sobol, top_k = 1L,
                                            score = score,
                                            rb_scoring = rb_scoring)
  base_params$mus[theta_names] <- theta_unc_hat
  base_params
}

# Random-search MVT mixture initialization for AMIS.
# Finds the top-K Sobol candidates and uses them as the K component means,
# giving genuinely diverse coverage rather than jittered copies of one point.
# The shared Sigma is taken from the heuristic joint init (same as AIS path).
init_mixture_params_random_search <- function(y, model_components, theta_names,
                                              eta_df, phi_min, phi_max,
                                              n_sobol, K) {
  base_params      <- init_joint_params(y, model_components, theta_names, eta_df)
  top_k_unc        <- find_theta_random_search(y, model_components,
                                               phi_min, phi_max,
                                               n_sobol = n_sobol, top_k = K)
  # top_k_unc is K x d_theta; embed into joint (theta + eta) parameter space
  param_names <- names(base_params$mus)
  eta_names   <- setdiff(param_names, theta_names)

  weights    <- rep(1 / K, K)
  mus_list   <- vector("list", K)
  Sigma_list <- vector("list", K)

  for (k in seq_len(K)) {
    mu_k                  <- base_params$mus
    mu_k[theta_names]     <- top_k_unc[k, ]
    names(mu_k)           <- param_names
    mus_list[[k]]         <- mu_k
    Sigma_list[[k]]       <- base_params$Sigma
    rownames(Sigma_list[[k]]) <- colnames(Sigma_list[[k]]) <- param_names
  }

  list(
    K          = K,
    weights    = weights,
    mus_list   = mus_list,
    Sigma_list = Sigma_list,
    df         = base_params$df
  )
}

# ---------------------------------------------------------------------------
# Low-level MVT helpers (no mvtnorm dependency)
# ---------------------------------------------------------------------------

# Sample n rows from MVT(mu, Sigma, df) given upper Cholesky R (R'R = Sigma).
rmvt_chol <- function(n, mu, chol_R, df) {
  d     <- length(mu)
  eps   <- matrix(stats::rnorm(n * d), n, d)
  Z     <- eps %*% chol_R                       # ~ MVN(0, Sigma)
  chi2  <- stats::rchisq(n, df = df)
  devs  <- Z / sqrt(chi2 / df)                  # zero-centred MVT deviations
  sweep(devs, 2, mu, "+")
}

# Log-density of MVT(mu, Sigma, df) for N rows, given:
#   devs   : N x d matrix of (x - mu) deviations
#   chol_R : upper Cholesky of Sigma
#   df     : degrees of freedom
# Returns a length-N numeric vector.
ldmvt_chol <- function(devs, chol_R, df) {
  d           <- ncol(devs)
  z           <- forwardsolve(t(chol_R), t(devs))      # d x N
  mahal       <- colSums(z^2)
  log_det_R   <- sum(log(diag(chol_R)))                 # = 0.5 * log|Sigma|
  log_const   <- lgamma((df + d) / 2) - lgamma(df / 2) -
                 (d / 2) * log(df * pi) - log_det_R
  log_const - ((df + d) / 2) * log(1 + mahal / df)
}

# Draw from Joint Proposal
draw_from_joint_proposal <- function(N, proposal_params,
                                     theta_names, phi_min, phi_max,
                                     antithetic = TRUE,
                                     chol_Sigma = NULL) {

  df          <- proposal_params$df
  mu          <- proposal_params$mus
  param_names <- names(mu)
  d           <- length(mu)

  # Pre-compute Cholesky of Sigma if not supplied by the caller.
  # When supplied, this decomposition is re-used across calls within one
  # AIS iteration, eliminating the double Cholesky that mvtnorm::rmvt and
  # mvtnorm::dmvt would each compute independently.
  if (is.null(chol_Sigma)) chol_Sigma <- chol(proposal_params$Sigma)

  # Sample from MVT(mu, Sigma, df) using the Cholesky factor:
  #   z = eps %*% R,  eps ~ N(0,I)  =>  z ~ MVN(0, Sigma)
  #   x = mu + z / sqrt(chi2/df),  chi2 ~ chi^2(df)
  N_half    <- ceiling(N / 2)
  eps       <- matrix(stats::rnorm(N_half * d), N_half, d)
  Z_half    <- eps %*% chol_Sigma                      # N_half x d, ~ MVN(0, Sigma)
  chi2      <- stats::rchisq(N_half, df = df)
  devs_half <- Z_half / sqrt(chi2 / df)               # zero-centred MVT deviations

  if (antithetic) {
    devs <- rbind(devs_half, -devs_half)[1:N, , drop = FALSE]
  } else {
    devs <- devs_half[1:N, , drop = FALSE]
  }
  samps_unc <- sweep(devs, 2, mu, "+")
  colnames(samps_unc) <- param_names

  # Log-density via ldmvt_chol — one triangular solve, no re-decomposition.
  log_density_unc <- ldmvt_chol(devs, chol_Sigma, df)

  # Split into theta and eta columns
  theta_cols <- which(param_names %in% theta_names)
  eta_cols   <- which(!param_names %in% theta_names)
  theta_unc  <- samps_unc[, theta_cols, drop = FALSE]
  eta_free   <- samps_unc[, eta_cols,   drop = FALSE]

  # Transform theta to constrained space; adjust log-density by Jacobian
  trans_res   <- transform_unconstrained_to_theta(theta_unc, theta_names, phi_min, phi_max)
  theta_con   <- trans_res$theta
  log_density <- log_density_unc - trans_res$log_jac

  # Handle seasonal sum-to-zero constraint for eta
  s_cols_free <- grep("^s\\d+$", colnames(eta_free), value = TRUE)
  if (length(s_cols_free) > 0) {
    last_s  <- -rowSums(eta_free[, s_cols_free, drop = FALSE])
    eta_con <- cbind(eta_free, last_s)
    m_idx   <- max(as.integer(sub("s", "", s_cols_free))) + 1
    colnames(eta_con)[ncol(eta_con)] <- paste0("s", m_idx)
  } else {
    eta_con <- eta_free
  }

  list(
    theta       = theta_con,
    eta         = eta_con,
    theta_unc   = theta_unc,
    eta_free    = eta_free,
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
