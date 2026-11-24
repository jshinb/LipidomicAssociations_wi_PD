# Reproduce the 252-lipid, adjusted-BMI GBA1-effect forest plots.
# In RStudio/Positron: set the working directory to the repo root, then source.

module_dir <- file.path(getwd(), "scripts_forestplot", "forestplot_252lipids")
if (!dir.exists(module_dir)) {
  stop("Set the working directory to the repository root before sourcing this script.")
}
project_root <- normalizePath(getwd())

source(file.path(module_dir, "01_forestplot_config.R"), local = FALSE)
source(file.path(module_dir, "02_forestplot_helpers.R"), local = FALSE)
source(file.path(module_dir, "03_forestplot_data.R"), local = FALSE)
source(file.path(module_dir, "04_build_forestplot.R"), local = FALSE)

results <- load_forestplot_data(project_root)
intervals <- prepare_forestplot_intervals(results, is_or = FALSE)
forestplot <- assemble_forestplot(intervals)
forestplot_with_legend <- add_forestplot_legend(forestplot)

output_dir <- file.path(project_root, "results", "forestplots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
date_tag <- format(Sys.Date(), "%Y-%m-%d")
base_path <- file.path(output_dir, paste0("forestplot_252lipids_gba1eff_PD_", date_tag, ".pdf"))
legend_path <- file.path(output_dir, paste0("forestplot_252lipids_gba1eff_PD_wi_legend_", date_tag, ".pdf"))

ggplot2::ggsave(base_path, plot = forestplot, width = 33.8, height = 17.8, units = "cm")
ggplot2::ggsave(legend_path, plot = forestplot_with_legend, width = 33.8, height = 17.8, units = "cm")
message("Wrote forest plots:\n  ", base_path, "\n  ", legend_path)
