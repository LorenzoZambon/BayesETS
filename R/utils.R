# TRUE for a single integer >= 1
is_count <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1 && x == round(x)
}

# Split n draws among the models, proportionally to their weights (the rounding
# remainder goes to the model with most draws)
draws_per_model <- function(model_weights, n) {
  n_k <- round(n * model_weights)
  n_k[which.max(n_k)] <- n - sum(n_k) + max(n_k)
  n_k
}
