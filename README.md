# Interval-Valued Repeated-Measures ANOVA with Permutation Inference

This repository implements an uncertainty-aware repeated-measures ANOVA
for interval-valued longitudinal data.
Each observation is represented by a center and radius, allowing inference
to jointly capture changes in location and uncertainty over time.

The proposed method constructs an ANOVA-style sum-of-squares decomposition
based on a weighted quadratic distance and performs inference using a
within-subject permutation test that is finite-sample exact under the
repeated-measures null hypothesis.

The code accompanying this repository reproduces the simulation results
and figures reported in the associated manuscript.

---

## Repository structure

nterval-RM-ANOVA/
├── R/ # Core functions (simulation, test statistic, permutation)
│ ├── sim_data_rm_interval_updated.R
│ ├── observe_F_rm_dw.R
│ ├── perm_within_subject_cr.R
│ ├── estimate_power_cr_parallel.R
│ ├── estimate_power_tau_grid.R
│ ├── make_sigma.R
│ └── utils.R
├── scripts/
│ └── run_power_study_updated.R # Main script to reproduce simulations + figures
├── figures/
│ └── poisson/ # Example output figures (PNG)
├── README.md
└── Dependent_ANOVA.Rproj


---

## Main script

To reproduce the simulation study and figures, run **one script**:

```r
scripts/run_power_study_updated.R
