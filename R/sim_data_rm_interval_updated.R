# Simulator: repeated-measures interval data ------------------------------
# Generates long data: id, time, center, radius
#
# Model (conceptual):
#   (C_it, R_it) = (muC_t, muR_t) + b_i + e_it
# where:
#   b_i is subject random intercept for both components (2D)
#   e_it is within-subject noise for both components (2D), optionally AR(1) over time per component
# Dependence between C and R is controlled by rho via the 2x2 covariance.
#
# Time effects:
#   - center_effect: numeric vector length T (added to muC_t)
#   - radius_effect: numeric vector length T (added to muR_t)
#
# Radii are forced nonnegative using pmax(., radius_floor).
#
# NEW (Jan 2026): Allow NON-NORMAL center noise too:
#   - center_design controls the marginal distribution for center noise
#   - radius_design controls the marginal distribution for radius noise
#   - if center_design is NULL, it defaults to radius_design (so both follow the same family)

source(file.path("R","utils.R"))
source(file.path("R","make_sigma.R"))

# ---- helper: map latent Gaussian Z -> standardized (mean 0, var 1) noise ----
.z_to_std_noise <- function(Z, design,
                            gamma_shape = 4, gamma_rate = 2,
                            poisson_lambda = 6,
                            t_df = 5) {

  design <- match.arg(design, c("normal", "norm_unif", "gamma", "poisson", "t"))

  if (design == "normal") {
    # Already mean 0, var 1
    return(Z)
  }

  U <- pnorm(Z)  # Gaussian copula -> Uniform(0,1)

  if (design == "norm_unif") {
    # Uniform(-sqrt(3), +sqrt(3)) has mean 0, var 1
    return(qunif(U, min = -sqrt(3), max = sqrt(3)))
  }

  if (design == "gamma") {
    # Standardize Gamma(shape, rate) to mean 0, var 1
    G <- qgamma(U, shape = gamma_shape, rate = gamma_rate)
    mG <- gamma_shape / gamma_rate
    vG <- gamma_shape / (gamma_rate^2)
    return((G - mG) / sqrt(vG))
  }

  if (design == "poisson") {
    if (poisson_lambda <= 0) stop("poisson_lambda must be > 0.")
    K <- qpois(U, lambda = poisson_lambda)
    return((K - poisson_lambda) / sqrt(poisson_lambda))
  }

  # Student t (heavy tails), standardized to var 1
  if (t_df <= 2) stop("t_df must be > 2 for finite variance.")
  X <- qt(U, df = t_df)
  # Var(t_df) = df/(df-2). Scale to var 1:
  return(X / sqrt(t_df / (t_df - 2)))
}

sim_data_rm_interval <- function(n = 50,
                                 T = 4,
                                 rho_CR = 0.5,
                                 sd_subject = c(0.7, 0.4),     # SD of subject intercept for (C,R)
                                 sd_noise   = c(1.0, 0.6),     # SD of within-subject noise for (C,R)
                                 ar1 = 0.0,                    # AR(1) across time applied per component
                                 center_effect = NULL,
                                 radius_effect = NULL,
                                 radius_floor = 0,

                                 # NEW: choose marginal families separately
                                 radius_design = c("norm_unif","gamma","poisson"),
                                 center_design = NULL,         # NULL -> same as radius_design; or set to "normal"/"gamma"/...

                                 # shared shape/params for non-normal families
                                 gamma_shape = 4,
                                 gamma_rate  = 2,
                                 poisson_lambda = 6,
                                 t_df = 5) {

  radius_design <- match.arg(radius_design)

  # if center_design not provided, match radius_design (common use-case)
  if (is.null(center_design)) {
    # radius_design is one of: norm_unif/gamma/poisson; map to same label
    center_design <- radius_design
  } else {
    # allow explicit center_design including "normal" and "t"
    center_design <- match.arg(center_design, c("normal","norm_unif","gamma","poisson","t"))
  }

  if (T < 2) stop("T must be >= 2.")
  if (length(sd_subject) != 2 || length(sd_noise) != 2) stop("sd_subject and sd_noise must have length 2.")
  if (!is.finite(ar1) || abs(ar1) >= 1) stop("ar1 must be in (-1,1).")
  if (!is.finite(rho_CR) || abs(rho_CR) >= 1) stop("rho_CR must be in (-1,1).")

  time_levels <- paste0("t", seq_len(T))
  if (is.null(center_effect)) center_effect <- rep(0, T)
  if (is.null(radius_effect)) radius_effect <- rep(0, T)
  if (length(center_effect) != T || length(radius_effect) != T) stop("center_effect and radius_effect must have length T.")

  # ----- Within-subject dependence (Gaussian copula backbone) -----
  Sigma_CR <- make_Sigma_CR(rho_CR)

  # AR(1) correlation across time
  Rtime <- outer(seq_len(T), seq_len(T), function(a,b) ar1^abs(a-b))

  # Build a correlation matrix for the stacked vector (C_1..C_T, R_1..R_T)
  Corr_time_C  <- 1.0 * Rtime
  Corr_time_R  <- 1.0 * Rtime
  Corr_time_CR <- Sigma_CR[1,2] * Rtime

  Corr_e <- rbind(
    cbind(Corr_time_C,  Corr_time_CR),
    cbind(Corr_time_CR, Corr_time_R)
  )

  # Simulate latent Gaussian noise with Corr_e (n x 2T)
  Ze <- rmvnorm_chol(n, Corr_e)

  # Split latent pieces
  Zc <- Ze[, 1:T, drop = FALSE]
  Zr <- Ze[, (T+1):(2*T), drop = FALSE]

  # Marginal transforms (mean 0, var 1) before scaling by sd_noise
  epsC_std <- .z_to_std_noise(
    Zc, design = center_design,
    gamma_shape = gamma_shape, gamma_rate = gamma_rate,
    poisson_lambda = poisson_lambda, t_df = t_df
  )

  epsR_std <- .z_to_std_noise(
    Zr, design = radius_design,
    gamma_shape = gamma_shape, gamma_rate = gamma_rate,
    poisson_lambda = poisson_lambda, t_df = t_df
  )

  # Scale to requested within-subject SDs
  epsC <- sd_noise[1] * epsC_std
  epsR <- sd_noise[2] * epsR_std

  # ----- Subject random intercepts (Gaussian) -----
  # Allow same rho_CR correlation between (bC, bR)
  Sigma_b2 <- diag(sd_subject, 2, 2) %*% Sigma_CR %*% diag(sd_subject, 2, 2)
  b <- rmvnorm_chol(n, Sigma_b2)  # n x 2
  bC <- b[, 1]
  bR <- b[, 2]

  # ----- Add time effects + intercepts -----
  muC <- center_effect
  muR <- radius_effect

  Cmat <- matrix(0, nrow = n, ncol = T)
  Rmat <- matrix(0, nrow = n, ncol = T)

  for (t in seq_len(T)) {
    Cmat[, t] <- muC[t] + bC + epsC[, t]
    Rmat[, t] <- muR[t] + bR + epsR[, t]
  }

  # Radii must be nonnegative
  Rmat <- pmax(Rmat, radius_floor)

  data.frame(
    id = rep(seq_len(n), each = T),
    time = factor(rep(time_levels, times = n), levels = time_levels),
    center = as.numeric(t(Cmat)),
    radius = as.numeric(t(Rmat))
  )
}
