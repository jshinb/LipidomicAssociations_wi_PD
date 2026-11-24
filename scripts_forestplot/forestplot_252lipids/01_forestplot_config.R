# Configuration for the 252-lipid GBA1-effect forest plot.
# Paths are resolved from the repository root by 05_run_forestplot.R.

FORESTPLOT_CLASS_LEVELS <- c(
  "TAG", "DAG", "CE", "Chol", "PC", "PC-O", "LPC", "PE", "PE-O",
  "LPE", "LPE-O", "PI", "Cer", "HexCer", "SM"
)

FORESTPLOT_COHORTS <- c(
  "GBA1eff: PD cases any vs. no mutations",
  "GBA1eff: controls any vs. no mutations"
)

FORESTPLOT_COLORS <- c(
  "GBA1eff: PD cases any vs. no mutations" = "#00BA38",
  "GBA1eff: controls any vs. no mutations" = "#619CFF"
)

forestplot_input_paths <- function(project_root) {
  list(
    lipid_classes = file.path(project_root, "Class252lipids.txt"),
    pd_cases = file.path(
      project_root, "results",
      "assoc_lipid_GBA1status_any_Species-Class_adjBMI_2026-10-08.txt"
    ),
    controls = file.path(
      project_root, "results",
      "assoc_test_gba1eff_Species-Class_Controls_adjBMI_2026-10-15.txt"
    )
  )
}
