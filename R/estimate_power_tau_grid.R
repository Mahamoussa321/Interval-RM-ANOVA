# R/estimate_power_tau_grid.R
# Grid over (w, n, tau) for power curves vs heterogeneity tau
# IMPORTANT: uses estimate_power_cr_parallel.R which is wired to
# sim_data_rm_interval_updated.R (supports center_design).

source(file.path("R", "estimate_power_cr_parallel.R"))

estimate_power_tau_grid <- function(n_list = c(10,30,50,70),
                                    w_list = list(0.3,0.6,1,"auto"),
                                    tau_list = seq(1.2, 2.0, by = 0.2),
                                    T = 4,
                                    reps = 200,
                                    B_perm = 2000,
                                    alpha = 0.05,
                                    design = "poisson",
                                    design_args = list(poisson_lambda = 6),
                                    
                                    # effect_strength scales the time-effect as tau increases
                                    effect_strength_center = 1.2,
                                    effect_strength_radius = 0.8,
                                    
                                    mode = c("center","radius"),
                                    seed = 1,
                                    
                                    # NEW: allow parallel reps inside estimate_power_cr
                                    parallel_reps = TRUE) {
  
  mode <- match.arg(mode)
  out <- list()
  idx <- 0
  
  for (w in w_list) {
    for (n in n_list) {
      for (tau in tau_list) {
        idx <- idx + 1
        
        if (mode == "center") {
          altC <- function(T) seq(0, effect_strength_center*(tau-1), length.out = T)
          altR <- function(T) rep(0, T)
        } else {
          altC <- function(T) rep(0, T)
          altR <- function(T) seq(0, effect_strength_radius*(tau-1), length.out = T)
        }
        
        df_one <- estimate_power_cr(
          n_list = c(n),
          T = T,
          reps = reps,
          B_perm = B_perm,
          alpha = alpha,
          alt_center = altC,
          alt_radius = altR,
          w = w,
          design = design,
          design_args = design_args,
          seed = seed + 10000*idx,
          verbose = FALSE,
          parallel_reps = parallel_reps
        )
        
        df_one$tau <- tau
        df_one$w <- if (identical(w, "auto")) "auto" else as.character(w)
        out[[idx]] <- df_one
      }
    }
  }
  
  do.call(rbind, out)
}
