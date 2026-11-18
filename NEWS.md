# BayesETS (development version)

* The Hessian at the mode no longer uses points where the integrand is not finite: one-sided differences there + bounded eigenvalues.
* AIS restarts the mode search (at most twice) from a draw for which the log-posterior is larger than its value at the mode + 2, as it means that the search had stopped too early.
* If AIS does not reach `min_ess` within `N_iter_max` iterations, the model keeps its estimate instead of getting zero weight, with a warning.
* An error in the fit of one model no longer stops `bets()`: the model gets zero weight, with a warning.

# BayesETS 0.1.0

* Initial CRAN submission.
