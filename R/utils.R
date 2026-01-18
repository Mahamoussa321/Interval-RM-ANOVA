
# Utilities ---------------------------------------------------------------

stopifnot_cols <- function(dat, cols) {
  miss <- setdiff(cols, names(dat))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

as_factor_safe <- function(x) {
  if (is.factor(x)) return(x)
  factor(x)
}

# Simple multivariate normal generator using Cholesky (base R only)
rmvnorm_chol <- function(n, Sigma) {
  p <- nrow(Sigma)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Z %*% L
}

