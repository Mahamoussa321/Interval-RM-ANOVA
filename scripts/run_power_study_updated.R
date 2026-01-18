# figures_power_size_all_dists_PERM_only_omega_FINAL_RUN.R
# ============================================================
# Creates TWO facet plots (PERM only) for EACH distribution:
#   (A) Type I error (size) at tau = 1  -> bar + CI, faceted by ω
#   (B) Power vs tau for tau > 1        -> line, facet grid (ω rows, n cols)
#
# Runs 3 designs (RM simulator):
#   1) gamma
#   2) poisson
#   3) norm_unif  (Uniform radii) with Normal centers forced via center_design="normal"
#
# Key fix:
#   - DO NOT pass radius_design in design_args (avoids duplicate match error)
#   - Use design="norm_unif" (not "normal_unif") for RM simulator
#
# Stability upgrades (borrowed from your successful "pub workflow"):
#   - L'Ecuyer-CMRG RNG for parallel reproducibility
#   - Conservative workers (detectCores()-2)
#   - future.globals.maxSize bump
#   - CHECKPOINT: save/load df RDS per DIST (resume after interrupt)
#   - tryCatch around estimate_power_tau_grid() so one DIST failure doesn't kill all
# ============================================================

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(ggplot2)
  library(dplyr)
  library(grid)  # for unit()
})

# -----------------------------
# RNG + Parallel plan (stable)
# -----------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(1)  # will be overwritten by SEED below, but keep deterministic init

options(future.globals.maxSize = 4 * 1024^3)  # 4GB safety for future globals
future::plan(multisession, workers = max(1L, parallel::detectCores() - 2L))
message("Workers: ", future::nbrOfWorkers())

# -----------------------------
# Source your pipeline
# -----------------------------
source(file.path("R", "estimate_power_tau_grid.R"))

# -----------------------------
# Global settings
# -----------------------------
N_LIST   <- c(5, 10, 30, 50)
W_LIST   <- list(0.3, 0.6, 1, "auto")
TAU_LIST <- c(1, seq(1.2, 2.0, by = 0.2))  # include tau=1 for Type I error
T_USE <- 3
MODE  <- "radius"   # "center" or "radius"
ALPHA <- 0.05
SEED  <- 1

REPS  <- 1000
BPERM <- 1000

# ensure deterministic master seed (and L'Ecuyer streams for workers)
set.seed(SEED)

# -----------------------------
# Colors (Okabe–Ito)
# -----------------------------
COL_BOOT <- "#0072B2"  # blue

# -----------------------------
# ω facet label (parsed)
# -----------------------------
make_omega_lab <- function(w_chr) {
  ifelse(w_chr == "auto", 'omega==omega^"*"',
         paste0("omega==", w_chr))
}

# -----------------------------
# Build design spec (RM simulator compatible)
#   - design controls radius family
#   - design_args controls parameters + center_design (when needed)
# -----------------------------
get_design_spec <- function(DIST) {
  if (DIST == "gamma") {
    list(
      design = "gamma",
      design_args = list(
        gamma_shape   = 2,
        gamma_rate    = 2,
        center_design = "gamma"   # ok (or NULL)
      )
    )
  } else if (DIST == "poisson") {
    list(
      design = "poisson",
      design_args = list(
        poisson_lambda = 6,
        center_design  = "poisson"  # ok (or NULL)
      )
    )
  } else if (DIST == "norm_unif") {
    # ENFORCE: Normal centers + Uniform radii
    # - radii: design="norm_unif"
    # - centers: center_design="normal"
    list(
      design = "norm_unif",
      design_args = list(
        center_design = "normal"
      )
    )
  } else {
    stop("Unknown DIST: ", DIST)
  }
}

DIST_LIST <- c( "norm_unif", "gamma","poisson")

# ============================================================
# MAIN LOOP: run + save two plots per distribution
#   - Checkpoint df per DIST in figures/<DIST>/df_PERM_<MODE>_<DIST>.rds
# ============================================================
for (DIST in DIST_LIST) {
  
  spec <- get_design_spec(DIST)
  DESIGN      <- spec$design
  DESIGN_ARGS <- spec$design_args
  
  message("\n==============================")
  message("Running DIST = ", DIST, " | design = ", DESIGN)
  message("==============================\n")
  
  # -----------------------------
  # Output dirs + checkpoint file
  # -----------------------------
  out_dir <- file.path("figures", DIST)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  rds_df <- file.path(out_dir, paste0("df_PERM_", MODE, "_", DIST, ".rds"))
  
  # -----------------------------
  # Compute (or resume)
  # -----------------------------
  if (file.exists(rds_df)) {
    message("Checkpoint found. Loading df from:\n  ", rds_df)
    df <- readRDS(rds_df)
  } else {
    
    df <- tryCatch(
      {
        estimate_power_tau_grid(
          n_list        = N_LIST,
          w_list        = W_LIST,
          tau_list      = TAU_LIST,
          T             = T_USE,
          reps          = REPS,
          B_perm        = BPERM,
          alpha         = ALPHA,
          design        = DESIGN,
          design_args   = DESIGN_ARGS,
          mode          = MODE,
          seed          = SEED,
          parallel_reps = TRUE
        )
      },
      interrupt = function(e) {
        message("Interrupted during DIST=", DIST, ". No df saved for this DIST yet.")
        return(NULL)
      },
      error = function(e) {
        message("ERROR during DIST=", DIST, ":\n  ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(df)) next
    
    # Save raw df checkpoint immediately (so plotting cannot lose it)
    saveRDS(df, rds_df)
    message("Saved df checkpoint:\n  ", rds_df)
  }
  
  df <- subset(df, method == "Permutation")
  df$method <- "Permutation"
  
  
  # ---- Labeling (ω and ~ n) ----
  df$n_lab <- paste0("n ~ ", df$n)
  df$n_lab <- factor(df$n_lab, levels = paste0("n ~ ", N_LIST))
  
  
  
  df$w <- ifelse(df$w == "auto", "auto", as.character(df$w))
  df$omega_lab <- make_omega_lab(df$w)
  
  # ============================================================
  # (A) Type I error (size): tau = 1
  # ============================================================
  df_size <- df %>%
    filter(tau == 1) %>%
    mutate(
      se = sqrt(pmax(0, estimate * (1 - estimate) / REPS)),
      lo = pmax(0, estimate - 1.96 * se),
      hi = pmin(1, estimate + 1.96 * se)
    )
  
  p_size <- ggplot(df_size, aes(x = factor(n), y = estimate)) +
    geom_hline(yintercept = ALPHA, linetype = 2, linewidth = 0.7) +
    geom_col(width = 0.65, fill = COL_BOOT, alpha = 0.85) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6) +
    facet_wrap(~ omega_lab, ncol = 2, labeller = labeller(omega_lab = label_parsed)) +
    labs(
      x = "Group size (~ n)",
      y = "False Rejection Rate"
  
    ) +
    coord_cartesian(ylim = c(0, .2)) +
    theme_bw(base_size = 13) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x       = element_text(angle = 25, hjust = 1),
      panel.grid.minor  = element_blank(),
      legend.position   = "top",
      legend.title      = element_blank(),
      panel.border      = element_rect(color = "gray40", fill = NA, linewidth = 1.0),
      strip.background  = element_rect(fill = "gray95", color = "gray40", linewidth = 1.0),
      strip.text        = element_text(face = "bold", size = 16),
      panel.spacing     = unit(0.7, "lines"),
      plot.margin       = margin(14, 22, 14, 14)
    )
  
  # ============================================================
  # (B) Power vs tau: tau > 1
  # ============================================================
  df_power <- df %>% filter(tau > 1)
  
  p_power <- ggplot(df_power, aes(x = tau, y = estimate, group = 1)) +
    geom_hline(yintercept = ALPHA, linetype = 2, linewidth = 0.7) +
    geom_line(linewidth = 0.9, color = COL_BOOT) +
    geom_point(size = 2.1, color = COL_BOOT) +
    facet_grid(
      omega_lab ~ n_lab,
      labeller = labeller(omega_lab = label_parsed)
    ) +
    labs(
      x = expression(paste("Heterogeneity factor  ", tau)),
      y = "Correct Rejection Rate"
    )+
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 16) +
    theme(
      panel.spacing.x   = unit(0.6, "lines"),
      panel.spacing.y   = unit(0.6, "lines"),
      strip.text.x = element_text(
        size = 16, face = "bold",
        lineheight = 0.9,
        margin = margin(t = 16, b = 16, l = 10, r = 10)
      ),
      strip.text.y = element_text(
        size = 16, face = "bold",
        lineheight = 0.9,
        margin = margin(t = 14, b = 14, l = 5, r = 5)
      ),
      strip.background  = element_rect(fill = "gray95", color = "gray40", linewidth = 1.0),
      axis.title        = element_text(size = 16, face = "bold"),
      axis.text         = element_text(size = 16),
      panel.border      = element_rect(color = "gray40", fill = NA, linewidth = 1.0),
      panel.grid.major  = element_line(color = "gray92", linewidth = 0.9),
      panel.grid.minor  = element_line(color = "gray96", linewidth = 0.8),
      plot.margin       = margin(16, 26, 16, 16),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 13)
    )
  
  print(p_size)
  print(p_power)
  
  # -----------------------------
  # Save figures
  # -----------------------------
  f1 <- file.path(out_dir, paste0("type1error_facets_PERM_", MODE, "_", DIST, ".png"))
  f2 <- file.path(out_dir, paste0("power_tau_facets_PERM_", MODE, "_", DIST, ".png"))
  
  ggsave(f1, p_size,  width = 10, height = 7, dpi = 300)
  ggsave(f2, p_power, width = 12, height = 8, dpi = 300)
  
  message("Saved:\n  ", f1, "\n  ", f2)
}

future::plan(future::sequential)
message("\nAll done.\n")
