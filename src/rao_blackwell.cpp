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

    for (int t = 0; t < L; t++) {
      // Forecast coefficient: yhat_t = forecast_coeff' * eta + forecast_det
      vec fc = coeff_l;
      double fd = det_l;
      if (trend) {
        fc += ph * coeff_b;
        fd += ph * det_b;
      }
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

      // Error deterministic part for state update
      double err_det = y(t) - fd;

      // State updates (coefficients propagation)
      vec lp_coeff = coeff_l;
      double lp_det = det_l;
      if (trend) {
        lp_coeff += ph * coeff_b;
        lp_det   += ph * det_b;
      }

      // l_{t+1} = l_t + phi*b_t + alpha*e_t
      //         = (l_t + phi*b_t) + alpha*(y_t - forecast)
      // coeff: lp_coeff + alpha*(-fc) = lp_coeff - alpha*fc
      vec new_coeff_l = lp_coeff - al * fc;
      double new_det_l = lp_det + al * err_det;

      vec new_coeff_b(n_eta, fill::zeros);
      double new_det_b = 0.0;
      if (trend) {
        // b_{t+1} = phi*b_t + beta*e_t
        new_coeff_b = ph * coeff_b - be * fc;
        new_det_b = ph * det_b + be * err_det;
      }

      if (seas) {
        // s_{j,t+m} = s_{j,t} + gamma*e_t (only for j = t%m)
        vec new_cs = coeff_s.row(sj).t() - ga * fc;
        coeff_s.row(sj) = new_cs.t();
        det_s(sj) = det_s(sj) + ga * err_det;
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
                            int L) {
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
  double half_L_log_psi0_nu0 = 0.5 * L * log(psi0 / nu0);

  // Precompute V0inv * eta0
  vec V0inv_eta0 = V0inv * eta0;
  double eta0_V0inv_eta0 = dot(eta0, V0inv_eta0);

  // Posterior covariance/mean matrices needed for reconstruction
  cube Vn_cube(n_eta, n_eta, N);
  mat  mu_n_mat(n_eta, N);
  vec  posterior_scale(N);  // psi_n = psi0 + quad for sigma^2 posterior

  for (int i = 0; i < N; i++) {
    mat XtXi = XtX_cube.slice(i);
    vec Xtyi = Xty_mat.col(i);
    double ytyi = yty_vec(i);

    // Vn = (V0inv + XtX)^{-1}  using Cholesky for stability
    mat M = V0inv + XtXi;
    mat Vn = inv_sympd(M);

    // Posterior mean: mu_n = Vn * (V0inv*eta0 + Xty)
    vec rhs = V0inv_eta0 + Xtyi;
    vec mu_n = Vn * rhs;

    Vn_cube.slice(i) = Vn;
    mu_n_mat.col(i) = mu_n;

    // Log determinant: log|I + X*V0*X'| = log|V0^{-1} + X'X| + log|V0|
    //                                    = log|M| + log|V0|
    //                                    = log|M| - log|V0inv|
    double log_det_M = log_det_sympd(M);
    double log_det_IpXVXt = log_det_M + log_det_V0;

    // Quadratic form: r = y~ - X*eta0
    // r'*(I+X*V0*X')^{-1}*r = rtr - Xtr' * Vn * Xtr
    // rtr = yty - 2*eta0'*Xty + eta0'*XtX*eta0
    double rtr = ytyi - 2.0 * dot(eta0, Xtyi) + dot(eta0, XtXi * eta0);

    // Xtr = Xty - XtX*eta0
    vec Xtr = Xtyi - XtXi * eta0;
    double quad_woodbury = dot(Xtr, Vn * Xtr);
    double quad = rtr - quad_woodbury;

    // Ensure numerical stability
    if (quad < 0.0) quad = 0.0;

    // Store posterior scale for sigma^2 sampling
    posterior_scale(i) = psi0 + quad;

    // Marginal Student-t log-likelihood
    // log p(y|theta) = lgamma((nu0+L)/2) - lgamma(nu0/2)
    //                - L/2 * log(nu0*pi)
    //                - 0.5 * log|Sigma_y|
    //                - (nu0+L)/2 * log(1 + quad/psi0)
    // where log|Sigma_y| = L*log(psi0/nu0) + log|I+X*V0*X'|

    double log_det_Sigma_y = half_L_log_psi0_nu0 + log_det_IpXVXt;

    log_ml(i) = lgamma_half_nu_L - lgamma_half_nu
              - half_L_log_nu_pi
              - 0.5 * log_det_Sigma_y
              - half_nu_L * log(1.0 + quad / psi0);
  }

  return List::create(
    _["log_marginal_lik"] = log_ml,
    _["Vn"] = Vn_cube,
    _["mu_n"] = mu_n_mat,
    _["posterior_scale"] = posterior_scale
  );
}
