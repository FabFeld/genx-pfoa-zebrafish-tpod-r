#' Generate sensitivity summary table and plot from annotated BMD results
#'
#' This function builds a pathway-level sensitivity summary (count, BMD quartiles,
#' and minimum BMD) from an annotated BMD data.frame (one row per gene/item)
#' and produces a ggplot showing pathway medians and IQRs per substance.
#'
#' @param res_anno data.frame Annotated BMD table. Must contain numeric `BMD.zSD`,
#'   a grouping column (default `pathway_name`) and a substance column
#'   (default `Substance`). Typically produced by merging per-substance bmd tables
#'   with a pathway annotation (ENSEMBL -> pathway_name).
#' @param output_dir character Directory where CSV and PDF will be saved. Created
#'   if it does not exist. If NULL, files are not written.
#' @param group_col character Name of the grouping column in `res_anno` (default: "pathway_name").
#' @param substance_col character Name of the substance column in `res_anno` (default: "Substance").
#' @param min_items integer Minimum number of items per group to include (default: 1).
#' @param save_csv logical Save the summary table as CSV to `output_dir` (default: TRUE).
#' @param save_plot logical Save the ggplot as PDF to `output_dir` (default: TRUE).
#' @param table_filename character Filename for CSV (default: "sensitivity_plot_table.csv").
#' @param plot_filename character Filename for PDF (default: "sensitivity_plots.pdf").
#' @param width numeric PDF width (inches). If NULL a sensible default is used.
#' @param height numeric PDF height (inches). If NULL a sensible default is used.
#'
#' @return A list with elements `table` (tibble) and `plot` (ggplot object).
#' @export
#' @examples
#' # generate_sensitivity_plot_from_res_anno(res_anno, output_dir = "comparison_plots")
generate_sensitivity_plot <- function(
  res_anno,
  output_dir = NULL,
  group_col = "pathway_name",
  substance_col = "Substance",
  min_items = 2,
  save_csv = TRUE,
  save_plot = TRUE,
  table_filename = "sensitivity_plot_table.csv",
  plot_filename = "sensitivity_plots.pdf",
  width = 8,
  height = NULL
) {
  if (!is.data.frame(res_anno)) stop("res_anno must be a data.frame")
  if (!group_col %in% names(res_anno)) stop("group_col not found in res_anno: ", group_col)
  if (!substance_col %in% names(res_anno)) stop("substance_col not found in res_anno: ", substance_col)
  if (!"BMD.zSD" %in% names(res_anno)) stop("res_anno must contain a numeric 'BMD.zSD' column")

  # ensure numeric BMD
  res_anno <- res_anno %>% dplyr::filter(!is.na(.data$BMD.zSD)) %>% dplyr::mutate(BMD.zSD = as.numeric(.data$BMD.zSD))

  # build summary table: per pathway and substance
  sensitivity_table <- res_anno %>%
    dplyr::group_by_at(c(group_col, substance_col)) %>%
    dplyr::summarize(
      nb_of_items = dplyr::n(),
      firstquartile = as.numeric(stats::quantile(BMD.zSD, 0.25, na.rm = TRUE)),
      secondquartile = as.numeric(stats::median(BMD.zSD, na.rm = TRUE)),
      thirdquartile = as.numeric(stats::quantile(BMD.zSD, 0.75, na.rm = TRUE)),
      min_BMD = as.numeric(min(BMD.zSD, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()%>% 
    dplyr::filter(nb_of_items >= min_items)

  if (nrow(sensitivity_table) == 0) stop("No groups remaining after filtering by min_items=", min_items)

  # order pathways by smallest median across substances (sensitive on top)
  pathway_order <- sensitivity_table %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarize(order_min = min(secondquartile, na.rm = TRUE)) %>%
    dplyr::arrange(order_min) %>%
    dplyr::pull(.data[[group_col]])

  sensitivity_table[[group_col]] <- factor(sensitivity_table[[group_col]], levels = pathway_order)

  # build ggplot: medians (secondquartile) with horizontal IQR errorbars
 

  plt <- ggplot2::ggplot(
    sensitivity_table,
    ggplot2::aes(
      x = secondquartile,
      y = .data[[group_col]],
      color = .data[[substance_col]],
      size = nb_of_items
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = firstquartile, xmax = thirdquartile), height = 0.2, linewidth = 0.6, alpha = 0.8) +
    (if (requireNamespace("ggsci", quietly = TRUE)) ggsci::scale_color_npg(alpha = 0.8) else ggplot2::scale_color_viridis_d()) +
    ggplot2::scale_x_log10(labels = tryCatch(format_log10_ticks(), error = function(e) scales::trans_format("log10", scales::math_format(10^.x)))) +
    ggplot2::labs(
      title = "Sensitivity plot by pathway",
      x = "BMD medians and IQRs",
      y = "Pathway",
      color = "Substance",
      size = "Number of items"
    ) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "right") +
    theme_bw()

  # save outputs if requested
  if (!is.null(output_dir) && (save_csv || save_plot)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(output_dir) && save_csv) {
    utils::write.csv(sensitivity_table, file = file.path(output_dir, table_filename), row.names = FALSE)
  }
  if (!is.null(output_dir) && save_plot) {
  if (is.null(height)) height <- max(6, 0.25 * length(unique(sensitivity_table[[group_col]])))
    ggplot2::ggsave(filename = file.path(output_dir, plot_filename), plot = plt, width = width, height = height)
  }

  invisible(list(table = sensitivity_table, plot = plt))
}

