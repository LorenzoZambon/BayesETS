#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// Build sufficient statistics X'X, X'y~, y~'y~ for the additive ETS
// Rao-Blackwellization.  Each particle has its own theta; the design
// matrix X(theta) and deterministic vector c(theta) are constructed in a
// single recursive pass per particle, and the sufficient statistics are
// accumulated on the fly so the full L x n_eta matrix never needs to be
// stored simultaneously.
//
// eta = (l0, [b0], [s0, ..., s_{m-1}])
// y_t = c_t(theta) + X_t(theta) * eta + eps_t
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List build_design_and_c_batch(NumericVector yR,
                              bool trend,
                              bool seas,
                              bool damped,
                              int m,
                              NumericMatrix paramsR) {

  vec y(yR.begin(), yR.size(), false);
  mat params(paramsR.begin(), paramsR.nrow(), paramsR.ncol(), false);

  int L = y.n_elem;
  int N = params.n_rows;

  int n_eta = 1;
  if (trend) n_eta += 1;
  if (seas) n_eta += m;

  vec alpha_vec = params.col(0);
  vec beta_vec, phi_vec, gamma_vec;
  int col_idx = 1;
  if (trend) {
    beta_vec = params.col(col_idx++);
    if (damped) {
      phi_vec = params.col(col_idx++);
    } else {
      phi_vec = vec(N, fill::ones);
    }
  }
  if (seas) {
    gamma_vec = params.col(col_idx);
  }

  // Output sufficient statistics
  cube XtX(n_eta, n_eta, N, fill::zeros);
  mat  Xty(n_eta, N, fill::zeros);
  vec  yty(N, fill::zeros);

  int s_offset = 1 + (trend ? 1 : 0);

  for (int i = 0; i < N; i++) {
    double al = alpha_vec(i);
    double be = trend ? beta_vec(i) : 0.0;
    double ph = trend ? phi_vec(i)  : 1.0;
    double ga = seas  ? gamma_vec(i) : 0.0;

    // State coefficient vectors over eta and deterministic parts.
    // coeff_l(j) = d(l_t)/d(eta_j), det_l = deterministic part of l_t
    vec coeff_l(n_eta, fill::zeros);
    coeff_l(0) = 1.0;
    double det_l = 0.0;

    vec coeff_b(n_eta, fill::zeros);
    double det_b = 0.0;
    if (trend) coeff_b(1) = 1.0;

    mat coeff_s;
    vec det_s;
    if (seas) {
      coeff_s = mat(m, n_eta, fill::zeros);
      det_s = vec(m, fill::zeros);
      for (int k = 0; k < m; k++)
        coeff_s(k, s_offset + k) = 1.0;
    }

    // Accumulators for sufficient statistics for this particle
    mat XtXi(n_eta, n_eta, fill::zeros);
    vec Xtyi(n_eta, fill::zeros);
    double ytyi = 0.0;

    // Preallocate temporaries outside the t-loop to avoid N*L heap allocations
    vec fc(n_eta), lp_coeff(n_eta), new_coeff_l(n_eta), new_coeff_b(n_eta);

    for (int t = 0; t < L; t++) {
      // Level + trend part of forecast — also used directly for state updates.
      fc = coeff_l;
      double fd = det_l;
      if (trend) {
        fc += ph * coeff_b;
        fd += ph * det_b;
      }
      // lp_coeff = coeff_l + ph*coeff_b is now in fc (pre-seasonal).
      // Capture it once here; avoids a redundant recomputation further down.
      lp_coeff = fc;
      double lp_det = fd;

      // Add seasonal contribution to complete the full forecast vector.
      int sj = 0;
      if (seas) {
        sj = t % m;
        fc += coeff_s.row(sj).t();
        fd += det_s(sj);
      }

      // y_tilde_t = y_t - c_t(theta) = y_t - forecast_det
      double yt = y(t) - fd;

      // Accumulate sufficient statistics: rank-1 updates
      XtXi += fc * fc.t();
      Xtyi += fc * yt;
      ytyi += yt * yt;

      // l_{t+1} = l_t + phi*b_t + alpha*e_t
      //         = (l_t + phi*b_t) + alpha*(y_t - forecast)
      // coeff: lp_coeff + alpha*(-fc) = lp_coeff - alpha*fc
      new_coeff_l = lp_coeff - al * fc;
      double new_det_l = lp_det + al * yt;

      double new_det_b = 0.0;
      if (trend) {
        // b_{t+1} = phi*b_t + beta*e_t
        new_coeff_b = ph * coeff_b - be * fc;
        new_det_b = ph * det_b + be * yt;
      }

      if (seas) {
        // s_{j,t+m} = s_{j,t} + gamma*e_t (only for j = t%m)
        coeff_s.row(sj) -= ga * fc.t();
        det_s(sj) = det_s(sj) + ga * yt;
      }

      coeff_l = new_coeff_l;
      det_l = new_det_l;
      if (trend) {
        coeff_b = new_coeff_b;
        det_b = new_det_b;
      }
    }

    XtX.slice(i) = XtXi;
    Xty.col(i) = Xtyi;
    yty(i) = ytyi;

  }

  return List::create(
    _["XtX"] = XtX,
    _["Xty"] = Xty,
    _["yty"] = yty,
    _["n_eta"] = n_eta,
    _["L"] = L
  );
}


// ---------------------------------------------------------------------------
// Marginal log-likelihood after analytically integrating out eta and sigma^2.
//
// Uses Woodbury identity / matrix determinant lemma to avoid L x L ops.
//
// y~ | theta ~ MVT(nu0, X*eta0, (psi0/nu0)*(I + X*V0*X'))
//
// log|I + X*V0*X'| = log|V0^{-1} + X'X| + log|V0|
// (I + X*V0*X')^{-1} = I - X*(V0^{-1} + X'X)^{-1}*X'
//
// Also returns posterior parameters Vn, mu_n, posterior_scale for
// posterior reconstruction of eta and sigma^2.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List marginal_likelihood_rb(arma::cube XtX_cube,
                            arma::mat  Xty_mat,
                            arma::vec  yty_vec,
                            arma::vec  eta0,
                            arma::mat  V0,
                            double nu0,
                            double psi0,
                            int L,
                            bool return_posterior = false) {
  int N = yty_vec.n_elem;
  int n_eta = eta0.n_elem;

  vec log_ml(N);

  // Precompute constants
  mat V0inv = inv_sympd(V0);
  double log_det_V0inv = log_det_sympd(V0inv);   // = -log|V0|
  double log_det_V0 = -log_det_V0inv;

  double half_nu_L = 0.5 * (nu0 + L);
  double half_nu   = 0.5 * nu0;
  double lgamma_half_nu_L = lgamma(half_nu_L);
  double lgamma_half_nu   = lgamma(half_nu);
  double half_L_log_nu_pi = 0.5 * L * log(nu0 * datum::pi);
  double L_log_psi0_nu0 = L * log(psi0 / nu0);

  // Precompute V0inv * eta0
  vec V0inv_eta0 = V0inv * eta0;

  // Only allocate posterior storage when the caller needs it
  cube Rn_cube;  // L_M^{-1} per particle (lower triangular); t(Rn)*Rn = Vn = M^{-1}
  mat  mu_n_mat;
  vec  posterior_scale;
  if (return_posterior) {
    Rn_cube.set_size(n_eta, n_eta, N);
    mu_n_mat.set_size(n_eta, N);
    posterior_scale.set_size(N);
  }

  for (int i = 0; i < N; i++) {
    mat XtXi = XtX_cube.slice(i);
    vec Xtyi = Xty_mat.col(i);
    double ytyi = yty_vec(i);

    // M = V0inv + XtX  (posterior precision)
    mat M = V0inv + XtXi;

    // Single lower Cholesky of M — yields log|M|, quadratic form, and
    // (on the posterior path) posterior mean and MVN sampling factor.
    mat L_M;
    if (!arma::chol(L_M, M, "lower")) {
      // M not positive definite: degenerate particle
      log_ml(i) = -datum::inf;
      if (return_posterior) posterior_scale(i) = datum::inf;
      continue;
    }

    // log|I + X*V0*X'| = log|M| + log|V0|;  log|M| = 2 * sum(log diag(L_M))
    double log_det_M = 2.0 * arma::sum(arma::log(L_M.diag()));
    double log_det_IpXVXt = log_det_M + log_det_V0;

    // Quadratic form: Xtr' * M^{-1} * Xtr = ||L_M^{-1} * Xtr||^2
    double rtr = ytyi - 2.0 * dot(eta0, Xtyi) + dot(eta0, XtXi * eta0);
    vec Xtr = Xtyi - XtXi * eta0;
    vec Lm_inv_Xtr = arma::solve(arma::trimatl(L_M), Xtr);
    double quad_woodbury = arma::dot(Lm_inv_Xtr, Lm_inv_Xtr);

    if (return_posterior) {
      // mu_n = M^{-1} * (V0inv*eta0 + Xty) via two triangular solves
      vec rhs = V0inv_eta0 + Xtyi;
      vec tmp = arma::solve(arma::trimatl(L_M), rhs);
      vec mu_n = arma::solve(arma::trimatu(L_M.t()), tmp);
      mu_n_mat.col(i) = mu_n;
      // Rn = L_M^{-1} (lower triangular). In R: crossprod(Rn, z) = t(Rn)*z
      // gives draws with covariance t(Rn)*Rn = L_M^{-T}*L_M^{-1} = M^{-1} = Vn.
      Rn_cube.slice(i) = arma::inv(arma::trimatl(L_M));
    }

    // quad >= 0 in exact arithmetic, but it is the difference of two terms
    // that become huge when the recursion explodes (theta outside the stable
    // region, long series).  If nearly all digits cancel, the value is noise:
    // mark the particle as invalid instead of clamping, so that rounding can
    // never produce a spuriously high likelihood.
    double quad = rtr - quad_woodbury;
    if (!std::isfinite(quad) || quad <= 1e-8 * rtr) {
      log_ml(i) = -datum::inf;
      if (return_posterior) posterior_scale(i) = datum::inf;
      continue;
    }

    // Marginal Student-t log-likelihood
    // log p(y|theta) = lgamma((nu0+L)/2) - lgamma(nu0/2)
    //                - L/2 * log(nu0*pi)
    //                - 0.5 * log|Sigma_y|
    //                - (nu0+L)/2 * log(1 + quad/psi0)
    // where log|Sigma_y| = L*log(psi0/nu0) + log|I+X*V0*X'|
    double log_det_Sigma_y = L_log_psi0_nu0 + log_det_IpXVXt;
    log_ml(i) = lgamma_half_nu_L - lgamma_half_nu
              - half_L_log_nu_pi
              - 0.5 * log_det_Sigma_y
              - half_nu_L * log(1.0 + quad / psi0);

    if (return_posterior) {
      posterior_scale(i) = psi0 + quad;
    }
  }

  if (return_posterior) {
    return List::create(
      _["log_marginal_lik"] = log_ml,
      _["Rn"]               = Rn_cube,
      _["mu_n"]             = mu_n_mat,
      _["posterior_scale"]  = posterior_scale
    );
  }
  return List::create(_["log_marginal_lik"] = log_ml);
}
