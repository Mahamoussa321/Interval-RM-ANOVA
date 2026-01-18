# Interval-Valued Repeated-Measures ANOVA with Permutation Inference

This repository implements an uncertainty-aware repeated-measures ANOVA for
**interval-valued longitudinal data**. Each observation is represented by a
**center** and a **radius**, allowing inference to jointly capture changes in
both location and uncertainty over time.

The proposed method constructs an ANOVA-style sum-of-squares decomposition
based on a **weighted quadratic distance** between intervals and performs
inference using a **within-subject permutation test**. Because the permutation
distribution is constructed under the exact null hypothesis, the resulting
p-values control Type I error in finite samples regardless of the marginal
distributions of the center and radius components.

The code in this repository reproduces the simulation results and figures
reported in the associated manuscript.

---

## Quick start (single script)

All simulations and figures can be reproduced by running **one script** from
the project root:

```r
source("scripts/run_power_study_updated.R")
Interval-RM-ANOVA/
├── R/                            # Core functions
│   ├── sim_data_rm_interval_updated.R
│   ├── observe_F_rm_dw.R
│   ├── perm_within_subject_cr.R
│   ├── estimate_power_cr_parallel.R
│   ├── estimate_power_tau_grid.R
│   ├── make_sigma.R
│   └── utils.R
├── scripts/
│   └── run_power_study_updated.R # Main script to reproduce simulations & figures
├── figures/
│   └── poisson/                  # Example output figures (PNG)
├── Dependent_ANOVA.Rproj
└── README.md

