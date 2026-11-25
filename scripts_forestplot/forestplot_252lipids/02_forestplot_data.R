# Read and assemble the adjusted association results used in this figure.

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

  read_cohort <- function(path, cohort_name, model_adj) {
    data.table::fread(path) |>
      dplyr::select(lipidID, beta, se, n, pval, fdr.p) |>
      dplyr::mutate(cohort = cohort_name) |>
      dplyr::mutate(adj = model_adj)
  }
  cohort_paths <- paths[FORESTPLOT_COHORTS]
  
  stopifnot(
    length(cohort_paths) == length(FORESTPLOT_COHORTS),
    length(cohort_paths) == length(FORESTPLOT_ADJ)
  )
  
  plot_data <- purrr::pmap_dfr(
    list(cohort_paths, FORESTPLOT_COHORTS, FORESTPLOT_ADJ),
    read_cohort
  ) |>
    dplyr::right_join(lipid_class, by = "lipidID") |>
    dplyr::relocate(lipidClass) |>
    dplyr::mutate(
      cohort = factor(cohort, levels = rev(FORESTPLOT_COHORTS))
    )
}

