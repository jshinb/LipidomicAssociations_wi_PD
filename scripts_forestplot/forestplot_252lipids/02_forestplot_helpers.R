# Shared helpers and forest-plot geoms.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggforestplot)
  library(ggpattern)
})

clean_lipid_id <- function(x) {
  x |>
    as.character() |>
    stringr::str_remove("\\*") |>
    stringr::str_remove_all(";\\d+") |>
    stringr::str_replace("_O-", "-O") |>
    stringr::str_replace("_", "") |>
    stringr::str_replace_all("_", "/")
}

generate_forest_col <- function(data, colors, is_or = FALSE,
                                 point_size = 2, stroke = 0.6,
                                 dodge_width = 0) {
  ggplot(data, aes(x = .estimate, y = lipidID)) +
    theme_forest() +
    geom_stripped_rows(odd = "grey85", even = "white", alpha = 0.2) +
    geom_vline(
      xintercept = if (is_or) 1 else 0,
      linetype = "solid", linewidth = 0.4, colour = "black"
    ) +
    geom_line(
      aes(x = beta, y = lipidID, colour = cohort, group = cohort),
      orientation = "y",
      position = position_dodge(width = dodge_width),
      linewidth = 0.5, alpha = 0.35
    ) +
    geom_linerange(
      aes(xmin = .xmin, xmax = .xmax, colour = cohort, alpha = p.signif, y = lipidID),
      position = position_dodge(width = dodge_width), linewidth = 0.5
    ) +
    scale_alpha_manual(values = c("p<0.05" = 0.9, "p>=0.05" = 0.2)) +
    geom_point(
      aes(x = beta, y = lipidID, colour = cohort, shape = fdr.signif, fill = fdr.signif),
      position = position_dodge(width = dodge_width),
      size = point_size, stroke = stroke
    ) +
    scale_color_manual(values = colors) +
    theme(axis.text.y = element_text(size = 5), legend.position = "none")
}
