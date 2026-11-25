# Modular 252-lipid forest plot

The numbered modules separate configuration, plotting helpers, input loading, plot construction, and the RStudio/Positron runner. The dated original script is unchanged.

## Run the real analysis

Set the RStudio/Positron working directory to the repository root and run:

```r
source("scripts_forestplot/forestplot_252lipids/05_run_forestplot.R")
```

The script reads `Class252lipids.txt` and the two adjusted association result files from `results/`, then writes PDFs to `results/forestplots/`.

## Smoke test with mock data

The synthetic fixtures contain invented lipid IDs and values; they are not study observations. From the repository root, run:

```r
source("tests/test_forestplot_252lipids.R")
```

The test creates temporary mock input files, checks loading, confidence intervals, and all five panels, then writes a temporary PDF and verifies it is non-empty. It does not modify the real analysis inputs or produce files in `results/`.

## Numbered modules

- `02_forestplot_helpers.R` — packages, lipid-ID normalization, and forest-plot geoms.
- `03_forestplot_data.R` — validates, reads, and joins input data.
- `04_build_forestplot.R` — confidence intervals, five lipid panels, combined plot, and legend.
- `05_run_forestplot.R` — config lipid-class order, cohorts, colors, and input paths then sources modules and writes the PDFs.

Required R packages: `dplyr`, `stringr`, `data.table`, `ggplot2`, `patchwork`, `ggforestplot`, and `ggpattern`.
