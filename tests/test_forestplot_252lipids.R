# Smoke test the modular forest-plot workflow using only synthetic fixture data.
# In RStudio/Positron: set the working directory to the repository root and source.

repo_root <- normalizePath(getwd())
module_dir <- file.path(repo_root, "scripts_forestplot", "forestplot_252lipids")
fixture_dir <- file.path(repo_root, "tests", "fixtures", "forestplot_252lipids")
stopifnot(dir.exists(module_dir), dir.exists(fixture_dir))

source(file.path(module_dir, "01_forestplot_config.R"))
source(file.path(module_dir, "02_forestplot_helpers.R"))
source(file.path(module_dir, "03_forestplot_data.R"))
source(file.path(module_dir, "04_build_forestplot.R"))

mock_root <- tempfile("forestplot-mock-")
dir.create(file.path(mock_root, "results"), recursive = TRUE)
mock_input_dir <- tempfile("forestplot-inputs-")
dir.create(mock_input_dir)
source(file.path(fixture_dir, "create_mock_inputs.R"), local = TRUE)

# The fixture generator writes into its evaluation directory; redirect/copy to
# the expected project-root locations for the normal data loader.
fixture_files <- c(
  "Class252lipids.txt",
  "assoc_lipid_GBA1status_any_Species-Class_adjBMI_2026-10-08.txt",
  "assoc_test_gba1eff_Species-Class_Controls_adjBMI_2026-10-15.txt"
)
# Regenerate in a temporary working directory so no repository files are changed.
old_wd <- setwd(mock_input_dir)
on.exit(setwd(old_wd), add = TRUE)
source(file.path(fixture_dir, "create_mock_inputs.R"), local = TRUE)
file.copy(file.path(mock_input_dir, fixture_files[[1]]), mock_root)
file.copy(file.path(mock_input_dir, fixture_files[[2]]), file.path(mock_root, "results"))
file.copy(file.path(mock_input_dir, fixture_files[[3]]), file.path(mock_root, "results"))
setwd(old_wd)

mock_results <- load_forestplot_data(mock_root)
stopifnot(nrow(mock_results) == 24L)
intervals <- prepare_forestplot_intervals(mock_results)
stopifnot(nrow(intervals) == 24L, all(intervals$.xmin <= intervals$.xmax))
panels <- build_forestplot_panels(intervals)
stopifnot(length(panels) == 5L, all(vapply(panels, nrow, integer(1)) > 0L))

plot <- assemble_forestplot(intervals)
legend_plot <- add_forestplot_legend(plot)
test_output <- file.path(mock_root, "smoke-test.pdf")
ggplot2::ggsave(test_output, plot = legend_plot, width = 25, height = 15, units = "cm")
stopifnot(file.exists(test_output), file.info(test_output)$size > 0)

message("Forest-plot smoke test passed. Mock output: ", test_output)
