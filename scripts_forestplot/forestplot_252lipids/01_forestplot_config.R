# Configuration for the 252-lipid GBA1-effect forest plot.
# Paths are resolved from the repository root by 05_run_forestplot.R.

FORESTPLOT_CLASS_LEVELS <- c(
  "TAG", "DAG", "CE", "Chol", "PC", "PC-O", "LPC", "PE", "PE-O",
  "LPE", "LPE-O", "PI", "Cer", "HexCer", "SM"
)

FORESTPLOT_COHORTS <- c(
  "GBA1: no or unknown mutations", 
  "PDeff: GBA1posPD vs. GBA1negHC"
)

FORESTPLOT_COLORS <- c(
  "GBA1: no or unknown mutations" = "#619CFF",
  "PDeff: GBA1posPD vs. GBA1negHC" = "#F8766D"
)

FORESTPLOT_ADJ = c('adjBMI','adjBMI')

forestplot_input_paths <- function(project_root) {
  list(
    lipid_classes = file.path(project_root, "Class252lipids.txt"),
    "GBA1: no or unknown mutations" = file.path(project_root,'results',"assoc_lipid_pd_Species-Class_GBA1neg_adjBMI_2026-10-08.txt"),
    "PDeff: GBA1posPD vs. GBA1negHC" = file.path(project_root,"results","assoc_test_PDposGBA_HCnegGBA_Species-Class_adjBMI_2026-11-24.txt")
  )
}
