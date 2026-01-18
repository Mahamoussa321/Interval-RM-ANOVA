
# Covariance builder for (center, radius) latent pair ---------------------
# Returns a 2x2 covariance for (C, R) with unit variances and correlation rho.
make_Sigma_CR <- function(rho = 0.5) {
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must be in (-1,1).")
  matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
}
