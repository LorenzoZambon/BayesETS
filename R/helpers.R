################################################################################
# HELPERS

#' Stacking weights
#'
#' @param log_lik_list List with one matrix of pointwise log-likelihoods
#'   (draws x time) per model.
#' @return Vector of stacking weights.
#' @keywords internal
compute_stacking_weights <- function(log_lik_list) {
  K <- length(log_lik_list)
  if (K == 0) stop("log_lik_list is empty")
  if (is.null(log_lik_list[[1]]) || !is.matrix(log_lik_list[[1]])) {
    stop("Each element of log_lik_list must be a matrix (samples x time)")
  }
  L <- ncol(log_lik_list[[1]])

  lpd_model <- matrix(nrow = L, ncol = K)
  for (k in seq_len(K)) {
    ll_mat <- log_lik_list[[k]]
    if (!is.matrix(ll_mat) || ncol(ll_mat) != L) {
      stop("All log_lik_list elements must be matrices with same number of columns")
    }
    S <- nrow(ll_mat)
    if (S <= 0) stop("Each log-likelihood matrix must have at least one row")
    max_ll <- apply(ll_mat, 2, max)
    sum_exp <- colSums(exp(ll_mat - matrix(max_ll, nrow = S, ncol = L, byrow = TRUE)))
    lpd_model[, k] <- max_ll + log(sum_exp) - log(S)
  }

  exp_lpd <- exp(lpd_model)

  obj_fun <- function(par) {
    w <- exp(par) / sum(exp(par))
    mix_dens <- as.vector(exp_lpd %*% w)
    -sum(log(pmax(mix_dens, 1e-300)))
  }

  opt <- stats::optim(rep(0, K), obj_fun, method = "BFGS")
  w_final <- exp(opt$par) / sum(exp(opt$par))

  w_final[w_final < 1e-3] <- 0
  if (sum(w_final) <= 0 || any(!is.finite(w_final))) {
    rep(1 / K, K)
  } else {
    w_final / sum(w_final)
  }
}

#' Label of an ETS model (e.g. "AAdN")
#'
#' @param model_components Vector (error, trend, season, damped).
#' @return Character string.
#' @keywords internal
ets_label <- function(model_components) {
  d <- if (model_components[[4]] == "TRUE" && model_components[[2]] == "A") "d" else ""
  paste0(model_components[[1]], model_components[[2]], d, model_components[[3]])
}

#' Future trajectories of the model combination
#'
#' Each model contributes a number of trajectories proportional to its weight.
#'
#' @param bets_fit The `fit` element of a `bets` object.
#' @param h Forecast horizon.
#' @param n_traj Number of trajectories.
#' @return Matrix of trajectories (n_traj x h).
#' @keywords internal
simulate_future_trajectories <- function(bets_fit, h = 10, n_traj = 1000) {
  n_models <- length(bets_fit$results)
  model_weights <- bets_fit$model_weights

  traj_list <- vector("list", n_models)
  n_traj_list <- round(n_traj * model_weights)
  n_traj_list[which.max(n_traj_list)] <- n_traj - sum(n_traj_list) + max(n_traj_list)

  for (i in seq_along(bets_fit$results)) {
    if (n_traj_list[i] > 0) {
      idxs <- sample(nrow(bets_fit$results[[i]]$thetas), size = n_traj_list[i], replace = TRUE)
      res_i <- bets_fit$results[[i]]
      traj_list[[i]] <- ets_future_traj(
        model_components = res_i$model_components,
        states = res_i$states[idxs, , drop = FALSE],
        params = res_i$thetas[idxs, , drop = FALSE],
        sigma2s = res_i$sigma2s[idxs],
        h = h
      )
    }
  }

  do.call(rbind, traj_list)
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


