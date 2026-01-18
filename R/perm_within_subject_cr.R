# Within-subject permutation for RM null ----------------------------------
# Under H0 (no time effect), time labels are exchangeable within each subject.
# We permute time within each id, recompute observe_F_rm_dw(), and build:
#   - p-value = mean(F* >= Fobs)
#   - cutoff = quantile(F*, 1-alpha)
#
# This gives an exact (finite-sample) test under exchangeability.

source(file.path("R","utils.R"))
source(file.path("R","observe_F_rm_dw.R"))   # NEW: notes-aligned F

perm_within_subject_cr <- function(dat_long, B = 5000, alpha = 0.05, seed = NULL, w = 1) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot_cols(dat_long, c("id","time","center","radius"))
  dat_long$id   <- as_factor_safe(dat_long$id)
  dat_long$time <- as_factor_safe(dat_long$time)
  
  obs <- observe_F_rm_dw(dat_long, w = w)    # NEW
  if (!is.finite(obs$F)) {
    return(list(F_obs = obs$F, df1 = obs$df1, df2 = obs$df2,
                p = NA_real_, cutoff = NA_real_, F_perm = numeric(0),
                note = obs$note))
  }
  
  ids <- levels(dat_long$id)
  Fstar <- numeric(B)
  
  # Pre-split rows by id for speed
  split_idx <- split(seq_len(nrow(dat_long)), dat_long$id)
  
  for (b in seq_len(B)) {
    dat_b <- dat_long
    
    for (s in ids) {
      idx <- split_idx[[s]]
      dat_b$time[idx] <- sample(dat_b$time[idx], length(idx), replace = FALSE)
    }
    
    dat_b$time <- factor(dat_b$time, levels = levels(dat_long$time))
    Fstar[b] <- observe_F_rm_dw(dat_b, w = w)$F   # NEW
  }
  
  Fstar <- Fstar[is.finite(Fstar)]
  cutoff <- as.numeric(stats::quantile(Fstar, probs = 1 - alpha, type = 8,
                                       names = FALSE, na.rm = TRUE))
  #pval <- mean(Fstar >= obs$F)
  pval <- (1 + sum(Fstar >= obs$F)) / (length(Fstar) + 1)
  
  list(F_obs = obs$F, df1 = obs$df1, df2 = obs$df2,
       p = pval, cutoff = cutoff, F_perm = Fstar, note = NULL)
}
