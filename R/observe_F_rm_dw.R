# Notes-aligned RM-ANOVA F for interval data using d_w^2 --------------------
# d_w^2([X],[Y]) = (C_X - C_Y)^2 + w (R_X - R_Y)^2
#
# RM decomposition (balanced):
#   SS_total = sum_{i,t} d_w^2( X_it, grand )
#   SS_treat = n * sum_{t} d_w^2( mean_t, grand )
#   SS_subj  = T * sum_{i} d_w^2( mean_i, grand )
#   SS_error = SS_total - SS_treat - SS_subj
#
# F = (SS_treat/(T-1)) / (SS_error/((T-1)(n-1)))

source(file.path("R","utils.R"))

observe_F_rm_dw <- function(dat_long, w = 1) {
  stopifnot_cols(dat_long, c("id","time","center","radius"))
  dat_long$id   <- as_factor_safe(dat_long$id)
  dat_long$time <- as_factor_safe(dat_long$time)
  
  tt  <- levels(dat_long$time)
  Tt  <- length(tt)
  ids <- levels(dat_long$id)
  n   <- length(ids)
  
  # Require balanced: one row per (id,time)
  dat_long <- dat_long[order(dat_long$id, dat_long$time), ]
  tab <- table(dat_long$id, dat_long$time)
  if (any(tab == 0)) stop("Unbalanced / missing time points. Make balanced first.")
  if (any(tab > 1)) stop("Multiple rows per (id,time). Aggregate first.")
  
  # --- Allow data-driven weight: w = sd(C)/sd(R) -------------------------
  if (is.character(w) && length(w) == 1 && w == "auto") {
    sdC <- stats::sd(dat_long$center, na.rm = TRUE)
    sdR <- stats::sd(dat_long$radius, na.rm = TRUE)
    if (!is.finite(sdC) || !is.finite(sdR) || sdR <= 0) {
      w <- 1
    } else {
      w <- sdC / sdR
    }
  }
  
  # Ensure w is a single finite numeric
  w <- suppressWarnings(as.numeric(w))
  if (!(length(w) == 1 && is.finite(w) && w >= 0)) {
    stop("w must be a single nonnegative number or w = 'auto'.")
  }
  w_used <- w
  
  # grand means
  grandC <- mean(dat_long$center, na.rm = TRUE)
  grandR <- mean(dat_long$radius, na.rm = TRUE)
  
  # time means
  time_means <- aggregate(cbind(center, radius) ~ time, data = dat_long, FUN = mean)
  # subject means
  subj_means <- aggregate(cbind(center, radius) ~ id, data = dat_long, FUN = mean)
  
  # SS_total
  SS_total <- sum((dat_long$center - grandC)^2 + w_used * (dat_long$radius - grandR)^2, na.rm = TRUE)
  
  # SS_treat
  SS_treat <- n * sum((time_means$center - grandC)^2 + w_used * (time_means$radius - grandR)^2, na.rm = TRUE)
  
  # SS_subj
  SS_subj  <- Tt * sum((subj_means$center - grandC)^2 + w_used * (subj_means$radius - grandR)^2, na.rm = TRUE)
  
  # SS_error
  SS_error <- SS_total - SS_treat - SS_subj
  
  df1 <- Tt - 1
  df2 <- (Tt - 1) * (n - 1)
  
  if (!is.finite(SS_error) || SS_error <= 0 || !is.finite(SS_treat)) {
    return(list(
      F = NA_real_, df1 = df1, df2 = df2,
      n = n, T = Tt, w_used = w_used,
      SS_total = SS_total, SS_treat = SS_treat, SS_subj = SS_subj, SS_error = SS_error,
      note = "Nonpositive/non-finite SS_error (numerical/degenerate case)."
    ))
  }
  
  MS_treat <- SS_treat / df1
  MS_error <- SS_error / df2
  Fstat <- MS_treat / MS_error
  
  list(
    F = Fstat, df1 = df1, df2 = df2,
    n = n, T = Tt, w_used = w_used,
    SS_total = SS_total, SS_treat = SS_treat, SS_subj = SS_subj, SS_error = SS_error,
    note = NULL
  )
}
