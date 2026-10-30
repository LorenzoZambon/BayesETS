# BayesETS 0.1.0

* Initial CRAN submission.
* `predict()` computes the forecast distribution in closed form, as a mixture of
  Student-t distributions: the forecasts are deterministic given the fit, without
  Monte Carlo error. Sample paths are simulated only with `return_paths = TRUE`,
  and `return_mixture = TRUE` returns the mixture.
