#' Generate Density Plots for BMD.zSD Values
#'
#' This function generates density plots for BMD.zSD values from a list of results
#' and saves the plots as a PDF file. It processes each result in the input list,
#' extracts the BMD.zSD values, and creates both absolute and proportional density
#' plots. The plots are combined into a single PDF file.
#'
#' @param results_list A named list of results, where each element contains data
#'   to extract BMD.zSD values. The names of the list are used for labeling.
#' @param output_dir A string specifying the directory where the output PDF file
#'   will be saved. The directory will be created if it does not exist.
#'
#' @details
#' The function performs the following steps:
#' - Creates the output directory if it does not exist.
#' - Iterates over the `results_list` to extract BMD.zSD values using the
#'   `extract_bmd_df` function.
#' - Filters out non-finite or missing BMD.zSD values.
#' - Generates density plots using the `dens_plot_multi` function for both
#'   absolute and proportional scales.
#' - Combines the plots into a single PDF file using `gridExtra::marrangeGrob`.
#'
#' @note
#' - If the `ggsci` package is available, it is used to generate a color palette
#'   for the plots. Otherwise, all plots will use black lines.
#' - If no valid BMD.zSD values are found in the input list, the function will
#'   stop with an error.
#'
#' @return
#' The function does not return a value. It saves the generated plots as a PDF
#' file in the specified output directory.
#'
#' @examples
#' \dontrun{
#' results <- list(
#'   "Sample1" = result1,
#'   "Sample2" = result2
#' )
#' generate_density_plots(results, "output/density_plots")
#' }
#'
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra marrangeGrob
#' @importFrom ggsci pal_npg
#' @export

generate_density_plots <- function(results_list, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  bmd_list <- list()
  for (nm in names(results_list)) {
    res <- results_list[[nm]]
    bmd_df <- extract_bmd_df(res)
    if (is.null(bmd_df)) { warning("No BMD table for '", nm, "'. Skipping."); next }
    if ("BMD.zSD" %in% names(bmd_df)) {
      bmd_df$BMD.zSD <- as.numeric(bmd_df$BMD.zSD)
      bmd_df <- bmd_df[is.finite(bmd_df$BMD.zSD) & !is.na(bmd_df$BMD.zSD), , drop = FALSE]
    } else { warning("Derived BMD.zSD not present for '", nm, "'. Skipping."); next }
    if (nrow(bmd_df) == 0) { warning("No finite BMD.zSD for '", nm, "'. Skipping."); next }
    bmd_list[[nm]] <- bmd_df
  }
  if (length(bmd_list) == 0) stop("No valid BMD tables available for density plotting.")
  colors <- if (requireNamespace("ggsci", quietly = TRUE)) ggsci::pal_npg()(length(bmd_list)) else rep("black", length(bmd_list))
  density_absolute <- dens_plot_multi(bmd_list, "BMD.zSD", show_elements = c("percentiles","gene_20","first_peak"), colors = colors, expand_factor = 0.3, proportional = FALSE, vline_size = 1, dens_line_width = 1)
  density_proportional <- dens_plot_multi(bmd_list, "BMD.zSD", show_elements = c("percentiles","gene_20","first_peak"), colors = colors, expand_factor = 0.3, proportional = TRUE, vline_size = 1, dens_line_width = 1)
  combined <- gridExtra::marrangeGrob(list(density_absolute, density_proportional), ncol = 1, nrow = 2)
  ggplot2::ggsave(filename = file.path(output_dir, "density_comparison.pdf"), plot = combined, width = 10, height = 9)
  message("Saved density comparison to ", file.path(output_dir, "density_comparison.pdf"))
}