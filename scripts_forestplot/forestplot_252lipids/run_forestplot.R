# Reproduce the 252-lipid, adjusted-BMI GBA1-effect forest plots.
# Run from any working directory with:
#   Rscript scripts_forestplot/forestplot_252lipids/run_forestplot.R

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
  module_dir <- dirname(script_path)
  project_root <- normalizePath(file.path(module_dir, "..", ".."))
} else {
  # Allows source() from the repository root in an interactive R session.
  project_root <- normalizePath(getwd())
  module_dir <- file.path(project_root, "scripts_forestplot", "forestplot_252lipids")
}

source(file.path(module_dir, "forestplot_config.R"), local = FALSE)
source(file.path(module_dir, "forestplot_helpers.R"), local = FALSE)
source(file.path(module_dir, "forestplot_data.R"), local = FALSE)
source(file.path(module_dir, "build_forestplot.R"), local = FALSE)

results <- load_forestplot_data(project_root)
intervals <- prepare_forestplot_intervals(results, is_or = FALSE)
forestplot <- assemble_forestplot(intervals)
forestplot_with_legend <- add_forestplot_legend(forestplot)

output_dir <- file.path(project_root, "results", "forestplots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
date_tag <- format(Sys.Date(), "%Y-%m-%d")

base_path <- file.path(
  output_dir, paste0("forestplot_252lipids_gba1eff_PD_", date_tag, ".pdf")
)
legend_path <- file.path(
  output_dir,
  paste0("forestplot_252lipids_gba1eff_PD_wi_legend_", date_tag, ".pdf")
)

ggplot2::ggsave(
  filename = base_path, plot = forestplot,
  width = 33.8, height = 17.8, units = "cm", device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  filename = legend_path, plot = forestplot_with_legend,
  width = 33.8, height = 17.8, units = "cm", device = grDevices::cairo_pdf
)

message("Wrote forest plots:\n  ", base_path, "\n  ", legend_path)
