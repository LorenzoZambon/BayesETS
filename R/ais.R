################################################################################
# ADAPTIVE IMPORTANCE SAMPLING (AIS)
#
# Samples the smoothing parameters only: the initial states and \sigma^2 are
# integrated analytically (conjugate prior). The proposal is a Student-t at the
# posterior mode, adapted while the ESS is below min_ess.

# AIS over the smoothing parameters, from z_start: log evidence and weighted
# draws of all iterations (theta, w, ml_res); failure: message, or NULL
integrate_ais <- function(log_g_fn, z_start, theta_names, ctrl,
                          log_g_mode = log_g_fn) {
  d <- length(theta_names)
  N_draw  <- resolve_by_d(ctrl$N_draw, d)
  min_ess <- if (is.null(ctrl$min_ess)) N_draw / 4 else resolve_by_d(ctrl$min_ess, d)

  ais <- adaptive_importance_sampling(log_g_fn, z_start,
                                      n_draw = N_draw, min_ess = min_ess,
                                      df = ctrl$is_df, scale = ctrl$is_scale,
                                      n_iter_max = ctrl$N_iter_max, lr = ctrl$lr,
                                      verbose = ctrl$verbose, log_g_mode = log_g_mode)
  proposal <- ais$proposal
  names(proposal$mus) <- theta_names
  rownames(proposal$Sigma) <- colnames(proposal$Sigma) <- theta_names

  list(
    log_evidence = ais$log_evidence,
    theta = do.call(rbind, lapply(ais$draw_evals, `[[`, "theta")),
    w = ais$w,
    ml_res = bind_ml_res(lapply(ais$draw_evals, `[[`, "ml_res")),
    ess = ais$ess,
    n_iter = ais$n_iter,
    proposal = proposal,
    timing = ais$timing,
    failure = if (ais$ess < min_ess) {
      sprintf("AIS: ESS = %.0f below min_ess = %.0f after %d iterations; the model gets zero weight",
              ais$ess, min_ess, ais$n_iter)
    }
  )
}


################################################################################
# AIS HELPERS

# AIS of g = exp(log g) on R^d (log_g_fn, z_start, log_g_mode as in
# adaptive_gh_quadrature()). Student-t proposal at the mode of log g, with scale
# matrix scale * H^{-1}; updated while ESS < min_ess, adding n_draw draws per
# iteration. The draws of all iterations are pooled and weighted with the
# mixture of the proposals.
adaptive_importance_sampling <- function(log_g_fn, z_start, n_draw, min_ess,
                                         df = 5, scale = 4,
                                         n_iter_max = 30, lr = 0.9,
                                         verbose = 0, log_g_mode = log_g_fn, ...) {
  lap <- laplace_mode(log_g_mode, z_start, ...)
  proposal <- list(mus = lap$zhat, Sigma = scale * tcrossprod(lap$Lmat), df = df)
  timing <- c(lap$timing, list(draws = 0, update = 0))

  proposals <- list()
  draw_evals <- list()
  Z <- NULL
  log_g <- NULL
  for (iter in seq_len(n_iter_max)) {
    # New draws, evaluated in one batch
    t0 <- proc.time()[3]
    Z_new <- draw_t_rqmc(n_draw, proposal$mus, proposal$Sigma, proposal$df)
    draw_eval <- log_g_fn(Z_new)
    timing$draws <- timing$draws + (proc.time()[3] - t0)
    proposals[[iter]]  <- c(proposal, list(n = nrow(Z_new)))
    draw_evals[[iter]] <- draw_eval
    Z <- rbind(Z, Z_new)
    log_g <- c(log_g, draw_eval$log_g)

    # Weights of the pooled draws
    log_w <- log_g - log_mix_density(Z, proposals)
    log_w[!is.finite(log_w)] <- -Inf
    lw_max <- max(log_w)
    if (!is.finite(lw_max)) {
      w <- NULL
      log_evidence <- -Inf
      ess <- 0
      break
    }
    w <- exp(log_w - lw_max)
    log_evidence <- lw_max + log(mean(w))
    w <- w / sum(w)
    ess <- 1 / sum(w^2)

    if (isTRUE(verbose >= 2)) {
      cat(sprintf("\n\nRao-Blackwellized AIS - iter %d\n", iter))
      cat(sprintf("\n ESS = %.1f\n", ess))
    }
    if (ess >= min_ess || iter == n_iter_max) break

    # Proposal update
    t0 <- proc.time()[3]
    proposal <- update_proposal(Z, w, proposal, lr = lr)
    timing$update <- timing$update + (proc.time()[3] - t0)
  }

  list(
    log_evidence = log_evidence,
    w = w,
    Z = Z,
    draw_evals = draw_evals,
    ess = ess,
    n_iter = iter,
    proposal = proposal,
    zhat = lap$zhat,
    timing = timing
  )
}

# n draws from a multivariate Student-t, from randomised Sobol points
draw_t_rqmc <- function(n, mu, Sigma, df) {
  d <- length(mu)
  u <- matrix(qrng::sobol(n, d = d + 1, randomize = "digital.shift"), ncol = d + 1)
  eps  <- stats::qnorm(u[, seq_len(d), drop = FALSE])
  chi2 <- stats::qchisq(u[, d + 1], df = df)
  devs <- (eps %*% chol(Sigma)) / sqrt(chi2 / df)
  # Drop draws from points exactly on 0 or 1
  devs <- devs[rowSums(!is.finite(devs)) == 0, , drop = FALSE]
  sweep(devs, 2, mu, "+")
}

# Log-density of the mixture of the proposals at the rows of Z (weights: share
# of the draws)
log_mix_density <- function(Z, proposals) {
  n_tot <- sum(vapply(proposals, `[[`, numeric(1), "n"))
  log_q <- vapply(proposals, function(p) {
    log(p$n / n_tot) + ldmvt_chol(sweep(Z, 2, p$mus, "-"), chol(p$Sigma), p$df)
  }, numeric(nrow(Z)))
  log_q <- matrix(log_q, nrow = nrow(Z))
  q_max <- apply(log_q, 1, max)
  q_max + log(rowSums(exp(log_q - q_max)))
}

# Binds the outputs of marginal_likelihood_rb() for several batches of particles
bind_ml_res <- function(ml_list) {
  if (length(ml_list) == 1L) return(ml_list[[1]])
  n_eta <- nrow(ml_list[[1]]$mu_n)
  n_tot <- sum(vapply(ml_list, function(r) length(r$posterior_scale), numeric(1)))
  list(
    log_marginal_lik = unlist(lapply(ml_list, function(r) as.numeric(r$log_marginal_lik))),
    Rn               = array(unlist(lapply(ml_list, `[[`, "Rn")), dim = c(n_eta, n_eta, n_tot)),
    mu_n             = do.call(cbind, lapply(ml_list, `[[`, "mu_n")),
    posterior_scale  = unlist(lapply(ml_list, function(r) as.numeric(r$posterior_scale))),
    final_coef       = array(unlist(lapply(ml_list, `[[`, "final_coef")), dim = c(n_eta, n_eta, n_tot)),
    final_const      = do.call(cbind, lapply(ml_list, `[[`, "final_const"))
  )
}

# Proposal update by weighted moment matching
update_proposal <- function(theta_unc, w, prev_params,
                            lr = 0.9, min_var = 1e-6,
                            lambda_shr = 0.1) {
  w <- w / sum(w)
  ess <- 1 / sum(w^2)

  new_mus <- colSums(w * theta_unc)

  centered <- sweep(theta_unc, 2, new_mus, "-")
  Sigma_new <- crossprod(centered * sqrt(w))

  diag_Sigma <- pmax(diag(Sigma_new), min_var)
  Sigma_new <- (1 - lambda_shr) * Sigma_new
  diag(Sigma_new) <- diag_Sigma

  lamb <- min(lr, ess / (100 + ess))

  out_mus <- lamb * new_mus + (1 - lamb) * prev_params$mus
  out_Sigma <- lamb * Sigma_new + (1 - lamb) * prev_params$Sigma

  list(mus = out_mus, Sigma = out_Sigma, df = prev_params$df)
}

# Log-density of a multivariate t, given the upper Cholesky factor of the scale
ldmvt_chol <- function(devs, chol_R, df) {
  d           <- ncol(devs)
  z           <- forwardsolve(t(chol_R), t(devs))
  mahal       <- colSums(z^2)
  log_det_R   <- sum(log(diag(chol_R)))
  log_const   <- lgamma((df + d) / 2) - lgamma(df / 2) - (d / 2) * log(df * pi) - log_det_R
  log_const - ((df + d) / 2) * log(1 + mahal / df)
}
