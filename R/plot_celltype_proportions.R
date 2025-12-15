#' Plot cell type proportions across visits
#'
#' @param object scLASER object
#' @param highlight_cell_types Optional character vector of cell type names
#'   to mark with an asterisk (e.g., interacting / changing types). Default:
#'   c("A", "B", "I", "J"). Set to NULL to disable highlighting.
#'
#' @return A ggplot object (invisibly)
#' @export
#'
plot_celltype_proportions <- function(object,
                                      highlight_cell_types = c("A", "B", "I", "J")) {
  stopifnot(inherits(object, "scLASER"))

  data <- object@metadata
  if (nrow(data) == 0) stop("metadata slot is empty in the scLASER object.")

  prop_data <- data %>%
    dplyr::group_by(subject_id, visit, disease, cell_type, sex) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(subject_id, visit, disease) %>%
    dplyr::mutate(
      n_total      = sum(count),
      prop         = count / n_total,
      proportion   = 100 * prop,
      sd_prop      = sqrt(prop * (1 - prop) / n_total),
      sd_percent   = 100 * sd_prop,
      # Works for any number of visits: V0, V1, V2, V3, ...
      visit_label  = paste0("V", visit),
      disease_label = ifelse(disease == 0, "Non-converter", "Converter"),
      sex_label     = ifelse(sex == 1, "Female", "Male")
    ) %>%
    dplyr::ungroup()

  # Mark selected cell types with an asterisk
  if (!is.null(highlight_cell_types)) {
    prop_data$cell_type2 <- ifelse(
      prop_data$cell_type %in% highlight_cell_types,
      paste0(prop_data$cell_type, "*"),
      prop_data$cell_type
    )
  } else {
    prop_data$cell_type2 <- prop_data$cell_type
  }

  p_prop <- ggplot2::ggplot(
    prop_data,
    ggplot2::aes(
      x     = visit_label,
      y     = proportion,
      group = subject_id,
      color = disease_label
    )
  ) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~cell_type2, nrow = 2) +
    ggplot2::scale_color_manual(
      values = c("Non-converter" = "#4DAF4A",
                 "Converter"     = "#984EA3")
    ) +
    ggplot2::labs(
      x     = "Visit",
      y     = "Proportion (%)",
      color = "Disease Status"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = ggplot2::element_text(size = 14),
      axis.text       = ggplot2::element_text(size = 12),
      axis.title      = ggplot2::element_text(size = 14),
      legend.text     = ggplot2::element_text(size = 12),
      legend.title    = ggplot2::element_text(size = 13)
    )

  print(p_prop)
  invisible(p_prop)
}
