# Transform results, assign lipids to five display panels, and assemble plots.

prepare_forestplot_intervals <- function(data, is_or = FALSE, plot.adj,
                                         FORESTPLOT_COHORTS,
                                         FORESTPLOT_CLASS_LEVELS,
                                         lipid_id_levels = NULL) {
  if(is.null(lipidID.levels)){
    lipid_id_levels <- data |>
      dplyr::distinct(lipidID2, lipidClass) |>
      dplyr::arrange(lipidClass, lipidID2) |>
      dplyr::pull(lipidID2) |>
      unique()
  }
  
  data |>
    dplyr::filter(adj == plot.adj) |>
    dplyr::mutate(
      lipidID = factor(lipidID2, levels = rev(lipid_id_levels)),
      lipidClass = factor(lipidClass, levels = FORESTPLOT_CLASS_LEVELS),
      cohort = factor(cohort, levels = rev(FORESTPLOT_COHORTS)),
      .xmin = if (is_or) exp(beta - 1.96 * se) else beta - 1.96 * se,
      .xmax = if (is_or) exp(beta + 1.96 * se) else beta + 1.96 * se,
      .estimate = if (is_or) exp(beta) else beta,
      fdr.signif = ifelse(fdr.p < 0.05, "p.fdr<0.05", "p.fdr>=0.05"),
      p.signif = ifelse(pval < 0.05, "p<0.05", "p>=0.05")
    ) |>
    dplyr::arrange(lipidClass, lipidID, cohort)
}

order_lipid_ids <- function(ids, lipid_id_levels) {
  ids <- unique(as.character(ids))
  ids[order(match(ids, rev(lipid_id_levels)), na.last = TRUE)]
}

select_panel_ids <- function(data, classes, carry_ids = character(),
                             limit = 50, lipid_id_levels) {
  class_ids <- data |>
    dplyr::filter(lipidClass %in% classes) |>
    dplyr::pull(lipidID) |>
    as.character() |>
    unique()
  ids <- order_lipid_ids(c(carry_ids, class_ids), lipid_id_levels)
  descending_ids <- rev(ids)
  selected <- head(descending_ids, limit)
  remaining <- if (length(selected) == length(descending_ids)) {
    character()
  } else {
    descending_ids[-seq_along(selected)]
  }
  list(selected = selected, remaining = remaining)
}

build_forestplot_panels <- function(intervals) {
  lipid_id_levels <- levels(intervals$lipidID)
  plot_data <- intervals |>
    dplyr::filter(cohort %in% FORESTPLOT_COHORTS)

  panel_1 <- plot_data |>
    dplyr::filter(lipidClass %in% c("TAG", "DAG")) |>
    dplyr::arrange(lipidClass, lipidID) |>
    dplyr::slice(rev(dplyr::row_number()))

  panel_2_ids <- select_panel_ids(plot_data, c("CE", "PC"), lipid_id_levels = rev(lipid_id_levels))
  panel_2 <- plot_data |> dplyr::filter(lipidID %in% panel_2_ids$selected)
  
  panel_3_ids <- select_panel_ids(
    plot_data, "PC-O", carry_ids = panel_2_ids$remaining,
    lipid_id_levels = rev(lipid_id_levels)
  )
  panel_3 <- plot_data |> dplyr::filter(lipidID %in% panel_3_ids$selected)
  
  panel_4_ids <- select_panel_ids(
    plot_data, c("LPC", "PE"), carry_ids = panel_3_ids$remaining,
    lipid_id_levels = rev(lipid_id_levels)
  )
  panel_4 <- plot_data |> dplyr::filter(lipidID %in% panel_4_ids$selected)
  panel_5_ids <- select_panel_ids(
    plot_data, c("PE-O", "LPE", "PI", "Cer", "SM"),
    carry_ids = panel_4_ids$remaining, limit = Inf,
    lipid_id_levels = rev(lipid_id_levels)
  )
  panel_5 <- plot_data |> dplyr::filter(lipidID %in% panel_5_ids$selected)

  list(panel_1, panel_2, panel_3, panel_4, panel_5)
}

assemble_forestplot <- function(intervals, colors = FORESTPLOT_COLORS,...) {
  panels <- build_forestplot_panels(intervals)
  plots <- lapply(panels, function(panel) generate_forest_col(panel, colors = colors))
  plot_limits <- range(intervals$.xmin, intervals$.xmax, na.rm = TRUE)

  Reduce(`+`, plots) +
    patchwork::plot_layout(ncol = 5, widths = rep(1, 5), axes = "collect", guides = "collect") &
    ggplot2::facet_grid(rows = ggplot2::vars(lipidClass), scales = "free_y", space = "free_y") &
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 6),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 6.125),
      strip.text.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 2.5, b = 5, l = 2.5)
    ) &
    ggplot2::scale_shape_manual(values = c("p.fdr<0.05" = 19, "p.fdr>=0.05" = 21)) &
    ggplot2::scale_fill_manual(values = c("p.fdr<0.05" = "black", "p.fdr>=0.05" = "white")) &
    ggplot2::coord_cartesian(xlim = plot_limits) &
    ggplot2::xlab("beta (SD units) +/- 95% CI") &
    ggplot2::ylab(NULL)
}

add_forestplot_legend <- function(plot) {
  plot &
    ggplot2::guides(
      color = ggplot2::guide_legend(reverse = TRUE, ncol = 1, override.aes = list(size = 1)),
      shape = ggplot2::guide_legend(reverse = FALSE, ncol = 1, override.aes = list(size = 1)),
      alpha = ggplot2::guide_legend(reverse = FALSE, ncol = 1, override.aes = list(size = 1))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.height = grid::unit(0.1, "cm"),
      legend.key.spacing.y = grid::unit(0.05, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    )
}
