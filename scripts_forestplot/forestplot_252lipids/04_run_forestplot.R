# Reproduce the 252-lipid, adjusted-BMI GBA1-effect forest plots.
# In RStudio/Positron: set the working directory to the repo root, then source.

#=============================================================================@
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

wd="~/Library/CloudStorage/OneDrive-Personal/00_Work/GermanLipidomics/"
setwd(wd)
#=============================================================================@

module_dir <- file.path(getwd(), "00_scripts", "forestplot_252lipids")
if (!dir.exists(module_dir)) {
  stop("Set the working directory to the repository root before sourcing this script.")
}
project_root <- normalizePath(getwd())

source(file.path(module_dir, "01_forestplot_config.R"), local = FALSE)
source(file.path(module_dir, "02_forestplot_helpers.R"), local = FALSE)
source(file.path(module_dir, "03_forestplot_data.R"), local = FALSE)
source(file.path(module_dir, "04_build_forestplot.R"), local = FALSE)

results <- load_forestplot_data(project_root)#assoc results
intervals <- prepare_forestplot_intervals(results, is_or = FALSE, plot.adj = 'adjBMI',
                                          FORESTPLOT_COHORTS = FORESTPLOT_COHORTS,
                                          FORESTPLOT_CLASS_LEVELS = FORESTPLOT_CLASS_LEVELS
                                          )
forestplot <- assemble_forestplot(intervals)
forestplot_with_legend <- add_forestplot_legend(forestplot)

output_dir <- file.path(project_root, "results", "forestplots_module")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
date_tag <- format(Sys.Date(), "%Y-%m-%d")
base_path <- file.path(output_dir, paste0("forestplot_252_PDeff_GBA1neg_GBA1_posPD_vs_negHC_", date_tag, ".pdf"))
legend_path <- file.path(output_dir, paste0("forestplot_252_PDeff_GBA1neg_GBA1_posPD_vs_negHC_", date_tag, ".pdf"))

ggplot2::ggsave(base_path, plot = forestplot, width = 33.8, height = 17.8, units = "cm")
ggplot2::ggsave(legend_path, plot = forestplot_with_legend, width = 33.8, height = 17.8, units = "cm")
message("Wrote forest plots:\n  ", base_path, "\n  ", legend_path)
