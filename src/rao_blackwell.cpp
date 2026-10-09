#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// Sufficient statistics X'X, X'y~ and y~'y~ of the additive ETS model
//   y_t = c_t(\theta) + X_t(\theta) \eta + e_t,  \eta = (l0, [b0], [m seasonal states]),
// with y~ = y - c, for each row (\theta) of params. X and c are built in one
// recursive pass, without storing X. With return_final, also the final states
// as affine functions of \eta: final_coef' \eta + final_const.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List build_design_and_c_batch(NumericVector yR,
                              bool trend,
                              bool seas,
                              bool damped,
                              int m,
                              NumericMatrix paramsR,
                              bool return_final = false) {

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

  // Final states: column j of final_coef.slice(i) holds the coefficients of state j
  cube final_coef;
  mat  final_const;
  if (return_final) {
    final_coef.zeros(n_eta, n_eta, N);
    final_const.zeros(n_eta, N);
  }

  int s_offset = 1 + (trend ? 1 : 0);

  for (int i = 0; i < N; i++) {
    double al = alpha_vec(i);
    double be = trend ? beta_vec(i) : 0.0;
    double ph = trend ? phi_vec(i)  : 1.0;
    double ga = seas  ? gamma_vec(i) : 0.0;

    // States as linear functions of eta: coefficients (coeff) and constant (det)
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

    // Sufficient statistics of this particle
    mat XtXi(n_eta, n_eta, fill::zeros);
    vec Xtyi(n_eta, fill::zeros);
    double ytyi = 0.0;

    // Temporaries, allocated once
    vec fc(n_eta), lp_coeff(n_eta), new_coeff_l(n_eta), new_coeff_b(n_eta);

    for (int t = 0; t < L; t++) {
      // Level + trend part of the forecast (also used in the state update)
      fc = coeff_l;
      double fd = det_l;
      if (trend) {
        fc += ph * coeff_b;
        fd += ph * det_b;
      }
      lp_coeff = fc;
      double lp_det = fd;

      // Seasonal part
      int sj = 0;
      if (seas) {
        sj = t % m;
        fc += coeff_s.row(sj).t();
        fd += det_s(sj);
      }

      // y~_t = y_t - c_t(\theta)
      double yt = y(t) - fd;

      // Rank-1 updates
      XtXi += fc * fc.t();
      Xtyi += fc * yt;
      ytyi += yt * yt;

      // l_{t+1} = l_t + \phi b_t + \alpha e_t, with e_t = y~_t - fc' \eta
      new_coeff_l = lp_coeff - al * fc;
      double new_det_l = lp_det + al * yt;

      double new_det_b = 0.0;
      if (trend) {
        // b_{t+1} = \phi b_t + \beta e_t
        new_coeff_b = ph * coeff_b - be * fc;
        new_det_b = ph * det_b + be * yt;
      }

      if (seas) {
        // s_{t+m} = s_t + \gamma e_t
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

    // Final states (l, [b,] seasonal states in order of use)
    if (return_final) {
      final_coef.slice(i).col(0) = coeff_l;
      final_const(0, i) = det_l;
      if (trend) {
        final_coef.slice(i).col(1) = coeff_b;
        final_const(1, i) = det_b;
      }
      if (seas) {
        for (int j = 0; j < m; j++) {
          int k = (L + j) % m;
          final_coef.slice(i).col(s_offset + j) = coeff_s.row(k).t();
          final_const(s_offset + j, i) = det_s(k);
        }
      }
    }
  }

  List out = List::create(
    _["XtX"] = XtX,
    _["Xty"] = Xty,
    _["yty"] = yty
  );
  if (return_final) {
    out["final_coef"] = final_coef;
    out["final_const"] = final_const;
  }
  return out;
}


// ---------------------------------------------------------------------------
// Final states of the posterior draws: final_coef_k' \eta_j + final_const_k,
// with k = idx(j) the particle (1-based) of draw j (column j of eta)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat final_states_rb(const arma::cube& final_coef,
                          const arma::mat& final_const,
                          const arma::mat& eta,
                          const arma::uvec& idx) {
  mat states(eta.n_rows, eta.n_cols);
  for (uword j = 0; j < eta.n_cols; j++) {
    uword k = idx(j) - 1;
    states.col(j) = final_coef.slice(k).t() * eta.col(j) + final_const.col(k);
  }
  return states;
}


// ---------------------------------------------------------------------------
// Log marginal likelihood p(y | \theta), with \eta and \sigma^2 integrated out:
//   y~ | \theta ~ t_{\nu_0}(X \eta_0, (\psi_0 / \nu_0) (I + X V_0 X'))
// computed with the Woodbury identity and the determinant lemma (no L x L
// matrices). Optionally returns the posterior parameters of \eta and \sigma^2.
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

  vec V0inv_eta0 = V0inv * eta0;

  // Posterior output, if requested
  cube Rn_cube;  // L_M^{-1}: t(Rn) Rn = M^{-1} = V_n
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

    // Cholesky of M: log|M|, quadratic form and posterior parameters
    mat L_M;
    if (!arma::chol(L_M, M, "lower")) {
      // M not positive definite: degenerate particle
      log_ml(i) = -datum::inf;
      if (return_posterior) posterior_scale(i) = datum::inf;
      continue;
    }

    // log|I + X V0 X'| = log|M| + log|V0|
    double log_det_M = 2.0 * arma::sum(arma::log(L_M.diag()));
    double log_det_IpXVXt = log_det_M + log_det_V0;

    // Xtr' M^{-1} Xtr = ||L_M^{-1} Xtr||^2
    double rtr = ytyi - 2.0 * dot(eta0, Xtyi) + dot(eta0, XtXi * eta0);
    vec Xtr = Xtyi - XtXi * eta0;
    vec Lm_inv_Xtr = arma::solve(arma::trimatl(L_M), Xtr);
    double quad_woodbury = arma::dot(Lm_inv_Xtr, Lm_inv_Xtr);

    if (return_posterior) {
      // mu_n = M^{-1} (V0^{-1} eta0 + X'y~)
      vec rhs = V0inv_eta0 + Xtyi;
      vec tmp = arma::solve(arma::trimatl(L_M), rhs);
      vec mu_n = arma::solve(arma::trimatu(L_M.t()), tmp);
      mu_n_mat.col(i) = mu_n;
      // Rn = L_M^{-1}: t(Rn) z has covariance V_n
      Rn_cube.slice(i) = arma::inv(arma::trimatl(L_M));
    }

    // quad >= 0, but it is a difference of two terms that are huge when the
    // recursion explodes (unstable \theta): if most digits cancel, the particle
    // is invalid (clamping at 0 would give a spurious high likelihood)
    double quad = rtr - quad_woodbury;
    if (!std::isfinite(quad) || quad <= 1e-8 * rtr) {
      log_ml(i) = -datum::inf;
      if (return_posterior) posterior_scale(i) = datum::inf;
      continue;
    }

    // Student-t log-density, with log|\Sigma_y| = L log(\psi_0 / \nu_0) + log|I + X V_0 X'|
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
