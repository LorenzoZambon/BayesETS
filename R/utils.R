# TRUE for a single integer >= 1
is_count <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1 && x == round(x)
}

# Indices (in their order) of the smallest set of largest weights holding at
# least 1 - mass of the total: the weights left out hold at most mass
heavy_set <- function(w, mass) {
  o <- order(w, decreasing = TRUE)
  n <- min(sum(cumsum(w[o]) < (1 - mass) * sum(w)) + 1L, length(w))
  sort(o[seq_len(n)])
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
