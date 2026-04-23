#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List RSS_vect_arma(NumericVector yR,
                   bool trend,
                   bool seas,
                   bool damped,
                   int m,
                   NumericMatrix init_statesR,
                   NumericMatrix paramsR,
                   bool return_residuals = false) {
  
  // Convert to Armadillo objects
  vec y(yR.begin(), yR.size(), false);
  mat init_states(init_statesR.begin(), init_statesR.nrow(), init_statesR.ncol(), false);
  mat params(paramsR.begin(), paramsR.nrow(), paramsR.ncol(), false);
  
  int L = y.n_elem;
  int N = params.n_rows;
  
  // Parameters
  vec alpha = params.col(0);
  vec beta, phi, gamma;
  
  if (trend) {
    beta = params.col(1);
    if (damped) {
      phi = params.col(2);
    } else {
      phi = vec(N, fill::ones);
    }
  }
  if (seas) {
    gamma = params.col(trend && damped ? 3 : trend ? 2 : 1);
  }
  
  // States
  mat l(N, L+1, fill::zeros);
  l.col(0) = init_states.col(0);
  
  mat b, s;
  if (trend) {
    b.set_size(N, L+1);
    b.col(0) = init_states.col(1);
  }
  if (seas) {
    s.set_size(N, L+m);
    for (int j=0; j<m; j++) {
      s.col(j) = init_states.col((trend ? 2 : 1) + m-1-j);;
    }
  }
  
  
  vec RSS(N, fill::zeros);
  mat residuals;
  if (return_residuals) {
    residuals.set_size(N, L);
    residuals.fill(0.0);
  }
  
  // Loop over time
  // Pre-allocate temporaries outside the loop to avoid N*L heap allocations
  vec ft(N);
  vec e(N);
  for (int t=0; t<L; t++) {
    ft = l.col(t);
    if (trend) ft += phi % b.col(t);  
    if (seas)  ft += s.col(t);
    
    // compute residuals and update RSS
    e = y[t] - ft;  
    if (return_residuals) {
      residuals.col(t) = e;
    }
    RSS += square(e);
    
    // update states
    l.col(t+1) = l.col(t) + alpha % e;
    if (trend) {
      l.col(t+1) += phi % b.col(t);
      b.col(t+1) = phi % b.col(t) + beta % e;
    }
    if (seas) {
      s.col(t+m) = s.col(t) + gamma % e;
    }
  }
  
  // Last states
  int ncols = 1 + (trend ? 1 : 0) + (seas ? m : 0);
  mat states(N, ncols, fill::zeros);
  states.col(0) = l.col(L);
  int col = 1;
  if (trend) {
    states.col(col) = b.col(L);
    col++;
  }
  if (seas) {
    for (int j=0; j<m; j++) {
      states.col(col) = s.col(L+j);
      col++;
    }
  }
  
  if (return_residuals) {
    return List::create(
      _["RSS"] = RSS,
      _["states"] = states,
      _["residuals"] = residuals
    );
  }
  return List::create(
    _["RSS"] = RSS,
    _["states"] = states
  );
}
