# Repeated-Measures ANOVA with Subject Bootstrap — Interval Data (Publication Project v2)

This project implements a subject-level (cluster) bootstrap for the *time* effect in repeated-measures ANOVA,
with simulation for **interval data** under three families: *Gamma*, *Poisson*, *Normal centers + Uniform radii*.
Scripts are numbered in the order you should run them.

## Order of running (scripts/)

1. `01_simulate_null_cutoff_per_family.R` — compute a 95% null cutoff for each family (Gamma/Poisson/Normal+Uniform).
2. `02_power_curve_per_family.R` — compute power curves per family (writes CSV + figures).
3. `03_three_families_example_analysis.R` — run one analysis per family with bootstrap p-values.
4. `04_sensitivity_wild_cluster.R` — optional: wild cluster bootstrap on one design (supplement).
5. `05_render_manuscript.R` — render the Quarto paper (PDF/HTML) using the results from steps 1–3.

## Quick start
- Run `install.R` to install required packages.
- Then run scripts in the order above.
- Figures go to `figures/`, results to `output/`, manuscript in `manuscript/`.
