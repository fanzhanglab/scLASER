#' Plot cell type proportions across visits
#'
#' @param object scLASER object
#' @param highlight_cell_types Optional character vector of cell type names
#'   to mark with an asterisk (e.g., interacting / changing types). Set to NULL
#'   to disable highlighting.
#' @param disease_var Column name in metadata encoding disease status
#' @param disease_labels Named character vector mapping disease values to labels
#'   (e.g., c("0" = "Non-responder", "1" = "Responder"))
#' @param disease_colors Named character vector mapping *labels* (not codes)
#'   to colors. Names must match the values produced by `disease_labels`.
#'
#' @return A ggplot object (invisibly)
#' @export
plot_celltype_proportions <- function(
    object,
    highlight_cell_types = NULL,
    disease_var = "disease",
    disease_labels,
    disease_colors
) {
  stopifnot(inherits(object, "scLASER"))

  data <- object@metadata
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("metadata slot is empty in the scLASER object.")
  }


  req <- c("subject_id", "visit", "cell_type", disease_var)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns in metadata: ", paste(miss, collapse = ", "))

  prop_data <- data %>%
    dplyr::group_by(
      .data$subject_id,
      .data$visit,
      .data[[disease_var]],
      .data$cell_type
    ) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data$subject_id, .data$visit, .data[[disease_var]]) %>%
    dplyr::mutate(
      n_total       = sum(.data$count),
      prop          = .data$count / .data$n_total,
      proportion    = 100 * .data$prop,
      sd_prop       = sqrt(.data$prop * (1 - .data$prop) / .data$n_total),
      sd_percent    = 100 * .data$sd_prop,
      visit_label   = paste0("V", .data$visit),
      disease_label = dplyr::recode(as.character(.data[[disease_var]]), !!!disease_labels)
    ) %>%
    dplyr::ungroup()

  # Mark selected cell types with an asterisk
  if (!is.null(highlight_cell_types)) {
    prop_data <- dplyr::mutate(
      prop_data,
      cell_type2 = dplyr::if_else(
        .data$cell_type %in% highlight_cell_types,
        paste0(.data$cell_type, "*"),
        .data$cell_type
      )
    )
  } else {
    prop_data$cell_type2 <- prop_data$cell_type
  }

  p_prop <- ggplot2::ggplot(
    prop_data,
    ggplot2::aes(
      x     = .data$visit_label,
      y     = .data$proportion,
      group = .data$subject_id,
      color = .data$disease_label
    )
  ) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(ggplot2::vars(.data$cell_type2), nrow = 2) +
    ggplot2::scale_color_manual(values = disease_colors) +
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
