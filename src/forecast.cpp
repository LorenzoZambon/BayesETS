#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ---------------------------------------------------------------------------
// Upper tail U(a) = P(T > a) and density of the standard t at a = 0, dz, ...,
// a_max, for cubic Hermite interpolation (error ~1e-10 with dz = 0.02)
// ---------------------------------------------------------------------------

struct TTable {
  double df, dz;
  int n;
  std::vector<double> U, f;
  TTable(double df_, double dz_ = 0.02, double a_max = 40) : df(df_), dz(dz_) {
    n = (int) std::round(a_max / dz) + 1;
    U.resize(n);
    f.resize(n);
    for (int i = 0; i < n; i++) {
      U[i] = R::pt(i * dz, df, 0, 0);
      f[i] = R::dt(i * dz, df, 0);
    }
  }
};

// CDF F and density d of the standard t at z: from the table (negative z by
// symmetry), exact beyond its range
inline void t_cdf_pdf(double z, const TTable& tb, double& F, double& d) {
  double a = std::fabs(z);
  double u = a / tb.dz;
  double U;
  if (u < tb.n - 1) {
    int i = (int) u;
    double t = u - i, t2 = t * t, t3 = t2 * t;
    double d0 = tb.f[i] * tb.dz, d1 = tb.f[i + 1] * tb.dz;
    U = (2 * t3 - 3 * t2 + 1) * tb.U[i] - (t3 - 2 * t2 + t) * d0 +
        (3 * t2 - 2 * t3) * tb.U[i + 1] - (t3 - t2) * d1;
    d = ((6 * t - 6 * t2) * (tb.U[i] - tb.U[i + 1]) + (3 * t2 - 4 * t + 1) * d0 +
         (3 * t2 - 2 * t) * d1) / tb.dz;
  } else {
    U = R::pt(a, tb.df, 0, 0);
    d = R::dt(a, tb.df, 0);
  }
  F = (z >= 0) ? 1 - U : U;
}


// ---------------------------------------------------------------------------
// t components of the forecast distribution of an additive ETS model at the
// particles idx (1-based): location f_j' m_T and scale sqrt(s2 (v_j + ||R_T f_j||^2))
// at horizons j = 1, ..., h (components x horizons), with m_T and R_T a column
// of final_mean and a slice of final_R, f_j' x = l + \phi_j b + (seasonal state
// of T + j), \phi_j = \phi + ... + \phi^j, and v_j = 1 + sum_{i < j} c_i^2,
// c_i = \alpha + \beta \phi_i + \gamma [i multiple of m]. States: l, [b,]
// [s1, ..., s_m], s1 used at T + 1. alpha, ..., s2: one value per particle in idx.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List ets_t_components(const arma::mat& final_mean,
                      const arma::cube& final_R,
                      const arma::uvec& idx,
                      const arma::vec& alpha,
                      const arma::vec& beta,
                      const arma::vec& phi,
                      const arma::vec& gamma,
                      const arma::vec& s2,
                      bool trend,
                      bool seas,
                      int h) {
  uword K = idx.n_elem;
  int s0 = trend ? 2 : 1;                                   // first seasonal state
  int m = seas ? (int) final_mean.n_rows - s0 : 1;
  mat location(K, h), scale(K, h);
  vec u(final_R.n_rows);

  for (uword k = 0; k < K; k++) {
    uword p = idx(k) - 1;
    const mat& R = final_R.slice(p);
    double phi_pow = 1, phi_j = 0, v = 1;
    for (int j = 0; j < h; j++) {                           // horizon j + 1
      phi_pow *= phi(k);
      phi_j += phi_pow;
      double loc = final_mean(0, p);
      u = R.col(0);                                         // R_T f_j
      if (trend) {
        loc += phi_j * final_mean(1, p);
        u += phi_j * R.col(1);
      }
      if (seas) {
        int sj = s0 + j % m;
        loc += final_mean(sj, p);
        u += R.col(sj);
      }
      location(k, j) = loc;
      scale(k, j) = std::sqrt(s2(k) * (v + dot(u, u)));
      double c = alpha(k) + beta(k) * phi_j + ((seas && (j + 1) % m == 0) ? gamma(k) : 0);
      v += c * c;
    }
  }
  return List::create(_["location"] = location, _["scale"] = scale);
}


// ---------------------------------------------------------------------------
// Quantiles probs of mixtures of location-scale t distributions with weights w
// and df degrees of freedom: one mixture per column of location and scale
// (components x mixtures); returns mixtures x probs. The quantile lies between
// the smallest and the largest quantile of the components: Newton steps on the
// mixture CDF inside this bracket, bisection when a step leaves it or shrinks
// too slowly (so convergence is ensured also for multimodal mixtures), until
// the step is below tol times the mean scale.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat t_mixture_quantiles(const arma::vec& w,
                              const arma::mat& location,
                              const arma::mat& scale,
                              double df,
                              const arma::vec& probs,
                              double tol = 1e-8) {
  uword K = location.n_rows, H = location.n_cols, P = probs.n_elem;
  if (w.n_elem != K || scale.n_rows != K || scale.n_cols != H)
    stop("w, location and scale must have matching dimensions");
  if (K == 0 || !(df > 0) || scale.min() <= 0)
    stop("at least one component, df > 0 and positive scales are needed");

  vec wn = w / accu(w);
  TTable tb(df);
  mat out(H, P);
  vec inv_s(K);

  for (uword h = 0; h < H; h++) {
    const double* loc = location.colptr(h);
    const double* sc = scale.colptr(h);
    double s_mean = 0;
    for (uword k = 0; k < K; k++) {
      inv_s(k) = 1 / sc[k];
      s_mean += wn(k) * sc[k];
    }
    double eps = tol * s_mean;

    for (uword j = 0; j < P; j++) {
      double p = probs(j);
      if (p <= 0 || p >= 1) {
        out(h, j) = (p <= 0) ? R_NegInf : R_PosInf;
        continue;
      }
      // Bracket and start: range and weighted mean of the component quantiles
      double qp = R::qt(p, df, 1, 0);
      double lo = R_PosInf, hi = R_NegInf, x = 0;
      for (uword k = 0; k < K; k++) {
        double q = loc[k] + sc[k] * qp;
        lo = std::min(lo, q);
        hi = std::max(hi, q);
        x += wn(k) * q;
      }
      if (hi - lo <= eps) {
        out(h, j) = x;
        continue;
      }

      double dx_old = hi - lo;
      for (int it = 0; it < 100; it++) {
        double F = 0, d = 0, Fk, dk;
        for (uword k = 0; k < K; k++) {
          t_cdf_pdf((x - loc[k]) * inv_s(k), tb, Fk, dk);
          F += wn(k) * Fk;
          d += wn(k) * dk * inv_s(k);
        }
        F -= p;
        if (F == 0) break;
        if (F < 0) lo = x; else hi = x;

        double xn = x - F / d;
        double dx = std::fabs(xn - x);
        if (!(xn > lo && xn < hi) || dx > 0.5 * dx_old) {
          xn = 0.5 * (lo + hi);
          dx = std::fabs(xn - x);
        }
        dx_old = dx;
        x = xn;
        if (dx <= eps || hi - lo <= eps) break;
      }
      out(h, j) = x;
    }
  }
  return out;
}
