# R/estimate_power_cr_parallel.R
# Power / Type-I simulation wrapper ----------------------------------------
# Notes-aligned: uses observe_F_rm_dw() (RM SS decomposition with d_w^2)
# PARALLEL VERSION: parallelize over Monte Carlo replicates (reps) using future.apply
#
# IMPORTANT:
#   This version uses sim_data_rm_interval_updated.R so CENTER noise can be non-normal.
#   We map the user-facing argument "design" -> simulator argument "radius_design".
#   To force center gamma too, pass in design_args = list(center_design="gamma", ...)

source(file.path("R","sim_data_rm_interval_updated.R"))
source(file.path("R","perm_within_subject_cr.R"))
source(file.path("R","observe_F_rm_dw.R"))

estimate_power_cr <- function(n_list = c(20,30,40,60,80),
                              T = 4,
                              reps = 1000,
                              B_perm = 2000,
                              alpha = 0.05,
                              alt_center = function(T) rep(0, T),
                              alt_radius = function(T) rep(0, T),
                              rho_CR = 0.5,
                              ar1 = 0.2,
                              w = 1,                       # omega weight
                              design = "norm_unif",         # we keep this name for your scripts
                              design_args = list(),         # forwarded to simulator
                              seed = 1,
                              verbose = TRUE,
                              parallel_reps = TRUE) {
  
  # Deterministic top-level seed (future.seed controls streams inside parallel lapply)
  set.seed(seed)
  
  if (parallel_reps) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel_reps=TRUE requires the 'future.apply' package. Install it or set parallel_reps=FALSE.")
    }
  }
  
  out_list <- vector("list", length(n_list))
  
  for (ii in seq_along(n_list)) {
    n <- n_list[ii]
    
    if (verbose) {
      message("Starting n=", n, " (reps=", reps, ", B_perm=", B_perm,
              ", w=", w, ", design=", design, ", T=", T, ")")
    }
    
    one_rep <- function(r_index) {
      
      ok <- FALSE
      tries <- 0
      
      while (!ok && tries < 10) {
        tries <- tries + 1
        
        # UPDATED simulator call:
        # - uses radius_design (NOT design)
        # - center_design is controlled via design_args; if NULL it defaults to radius_design
        dat <- do.call(sim_data_rm_interval, c(list(
          n = n,
          T = T,
          rho_CR = rho_CR,
          ar1 = ar1,
          center_effect = alt_center(T),
          radius_effect = alt_radius(T),
          radius_design = design
        ), design_args))
        
        obs <- tryCatch(observe_F_rm_dw(dat, w = w), error = function(e) NULL)
        
        if (is.null(obs) || is.null(obs$F) || is.null(obs$df1) || is.null(obs$df2)) {
          next
        }
        
        Fobs <- suppressWarnings(as.numeric(obs$F))
        df1  <- suppressWarnings(as.numeric(obs$df1))
        df2  <- suppressWarnings(as.numeric(obs$df2))
        
        ok <- (length(Fobs) == 1 && length(df1) == 1 && length(df2) == 1 &&
                 is.finite(Fobs) && is.finite(df1) && is.finite(df2) &&
                 df1 > 0 && df2 > 0)
      }
      
      if (!ok) {
        return(c(rej_perm = FALSE, rej_class = FALSE))
      }
      
      # Classic p-value (F reference)
      p_classic <- pf(Fobs, df1 = df1, df2 = df2, lower.tail = FALSE)
      
      # Permutation p-value (same statistic, RM-null)
      perm <- perm_within_subject_cr(dat, B = B_perm, alpha = alpha, w = w)
      p_perm <- perm$p
      
      c(rej_perm  = (p_perm < alpha),
        rej_class = (p_classic < alpha))
    }
    
    res_r <- if (parallel_reps) {
      future.apply::future_lapply(seq_len(reps), one_rep, future.seed = TRUE)
    } else {
      lapply(seq_len(reps), one_rep)
    }
    
    rej_perm  <- vapply(res_r, function(x) as.logical(x["rej_perm"]),  logical(1))
    rej_class <- vapply(res_r, function(x) as.logical(x["rej_class"]), logical(1))
    
    est_perm  <- mean(rej_perm)
    est_class <- mean(rej_class)
    se_perm   <- sqrt(est_perm  * (1 - est_perm)  / reps)
    se_class  <- sqrt(est_class * (1 - est_class) / reps)
    
    lo_perm   <- max(0, est_perm  - 1.96 * se_perm)
    hi_perm   <- min(1, est_perm  + 1.96 * se_perm)
    lo_class  <- max(0, est_class - 1.96 * se_class)
    hi_class  <- min(1, est_class + 1.96 * se_class)
    
    if (verbose) {
      message("n=", n,
              "  perm=", round(est_perm, 3),
              "  classic=", round(est_class, 3),
              "  (w=", w, ", design=", design, ")")
    }
    
    out_list[[ii]] <- data.frame(
      n = rep(n, 2),
      method = c("Permutation", "Classic F"),
      estimate = c(est_perm, est_class),
      se = c(se_perm, se_class),
      lo = c(lo_perm, lo_class),
      hi = c(hi_perm, hi_class),
      reps = reps,
      B_perm = B_perm,
      w = w,
      design = design
    )
  }
  
  do.call(rbind, out_list)
}
