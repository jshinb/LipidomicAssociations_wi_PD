# Read and assemble the adjusted association results used in this figure.
# The source script ultimately filters to adjBMI, so the unadjusted RData input
# and unrelated cohort result files are intentionally not loaded here.

load_forestplot_data <- function(project_root) {
  paths <- forestplot_input_paths(project_root)
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing)) {
    stop(
      "Required forest-plot input file(s) not found:\n  ",
      paste(missing, collapse = "\n  "),
      "\nSee the module README for expected input locations.",
      call. = FALSE
    )
  }

  lipid_class <- data.table::fread(paths$lipid_classes) |>
    dplyr::mutate(
      lipidClass = stringr::str_replace(lipidClass, "_O[-]", "-O"),
      lipidClass = factor(lipidClass, levels = FORESTPLOT_CLASS_LEVELS),
      lipidID2 = clean_lipid_id(lipidID)
    )

  read_cohort <- function(path, cohort_name) {
    data.table::fread(path) |>
      dplyr::select(lipidID, beta, se, n, pval, fdr.p) |>
      dplyr::mutate(cohort = cohort_name)
  }

  dplyr::bind_rows(
    read_cohort(paths$pd_cases, FORESTPLOT_COHORTS[[1]]),
    read_cohort(paths$controls, FORESTPLOT_COHORTS[[2]])
  ) |>
    dplyr::right_join(lipid_class, by = "lipidID") |>
    dplyr::relocate(lipidClass) |>
    dplyr::mutate(
      cohort = factor(cohort, levels = rev(FORESTPLOT_COHORTS)),
      adj = "adjBMI"
    )
}
