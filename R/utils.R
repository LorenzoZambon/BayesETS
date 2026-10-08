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

# Prints a character matrix with its column names as header, indented by two
# spaces: first column left-aligned, the others right-aligned
cat_table <- function(tab) {
  tab <- rbind(colnames(tab), tab)
  width <- apply(nchar(tab), 2, max)
  for (i in seq_len(nrow(tab))) {
    cells <- c(sprintf("%-*s", width[1], tab[i, 1]), sprintf("%*s", width[-1], tab[i, -1]))
    cat("  ", paste(cells, collapse = "  "), "\n", sep = "")
  }
}
