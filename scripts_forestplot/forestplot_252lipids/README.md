# Modular 252-lipid forest plot

This folder separates the original forest-plot workflow into configuration, shared plotting helpers, input loading, plot construction, and an executable driver. The existing dated script is left unchanged.

## Run

From the repository root:

```sh
Rscript scripts_forestplot/forestplot_252lipids/run_forestplot.R
```

The driver also resolves the repository root from its own script path, so it can be launched from another working directory when given the script path.

## Required inputs

Place these files at the listed paths relative to the repository root:

- `Class252lipids.txt`
- `results/assoc_lipid_GBA1status_any_Species-Class_adjBMI_2026-10-08.txt`
- `results/assoc_test_gba1eff_Species-Class_Controls_adjBMI_2026-10-15.txt`

The original script loads additional adjusted and unadjusted result files, but ultimately selects only `adjBMI` and the PD-cases-any-versus-no-mutations and controls-any-versus-no-mutations cohorts. This modular version loads only those inputs required for the displayed forest plots.

Expected association-result columns include `lipidID`, `beta`, `se`, `n`, `pval`, and `fdr.p`. Inputs are not included here; use the corresponding project-generated result files.

## Modules

- `forestplot_config.R` — lipid-class order, cohorts, colors, and input paths.
- `forestplot_helpers.R` — packages, lipid-ID normalization, and forest-plot geoms.
- `forestplot_data.R` — validates, reads, and joins the input data.
- `build_forestplot.R` — confidence intervals, five lipid panels, combined plot, and legend.
- `run_forestplot.R` — sources modules and writes the PDFs.

Outputs are written to `results/forestplots/` with the run date in each filename. The no-legend and legend PDFs are both produced. For reproducible comparisons, retain the input result files and record the R/package versions and run date; the displayed CI is beta +/- 1.96 * SE and significance thresholds are unchanged from the source script.

Required R packages: `dplyr`, `stringr`, `data.table`, `ggplot2`, `patchwork`, `ggforestplot`, and `ggpattern`.
