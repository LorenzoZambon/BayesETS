# integrate_ais() that always fails: with with_mocked_bindings(), the models
# fitted by AIS get zero weight
failing_ais <- function(...) {
  list(log_evidence = -Inf, failure = "AIS: mocked failure; the model gets zero weight")
}
