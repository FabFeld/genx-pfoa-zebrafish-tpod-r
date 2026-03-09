#' Build Sensitivity Summary Table
#'
#' This function generates a summary table from annotated BMD (Benchmark Dose) results.
#' It groups the data by specified columns, calculates summary statistics for the `BMD.zSD` column,
#' and filters groups based on the minimum number of items.
#'
#' @param res_anno A `data.frame` containing at least a `BMD.zSD` column and a grouping column.
#' @param group_col A `character` string specifying the column name to group by (default: `"pathway_name"`).
#' @param substance_col A `character` string specifying the column name for substances (default: `"Substance"`).
#' @param min_items An `integer` specifying the minimum number of items required in a group (default: `1`).
#'
#' @return A `data.frame` summarizing the sensitivity data with the following columns:
#' \itemize{
#'   \item `groupby`: The grouping column renamed from `group_col`.
#'   \item `Substance`: The substance column.
#'   \item `nb_of_items`: The number of items in each group.
#'   \item `firstquartile`: The first quartile of `BMD.zSD`.
#'   \item `secondquartile`: The median of `BMD.zSD`.
#'   \item `thirdquartile`: The third quartile of `BMD.zSD`.
#'   \item `min_BMD`: The minimum value of `BMD.zSD`.
#' }
#'
#' @details
#' The function filters out rows where `BMD.zSD` is `NA` and ensures that the `BMD.zSD` column is numeric.
#' It also checks for the existence of the specified grouping and substance columns in the input data.
#' Groups with fewer items than `min_items` are excluded from the output.
#'
#' @examples
#' # Example usage:
#' # build_sensitivity_table(res_anno, group_col = "pathway_name", substance_col = "Substance", min_items = 5)
#'
#' @importFrom dplyr filter mutate group_by_at summarize ungroup n
#' @importFrom stats quantile median
#' @export
## Build sensitivity summary table from annotated BMD results -----------------
build_sensitivity_table <- function(res_anno, group_col = "pathway_name", substance_col = "Substance", min_items = 1) {
  # res_anno: data.frame containing at least a BMD.zSD column and a pathway/group column
  if (!is.data.frame(res_anno)) stop("res_anno must be a data.frame")
  if (!"BMD.zSD" %in% names(res_anno)) stop("res_anno must contain a 'BMD.zSD' column")

  df <- res_anno %>% dplyr::filter(!is.na(BMD.zSD)) %>% dplyr::mutate(BMD.zSD = as.numeric(BMD.zSD))

  # ensure grouping columns exist
  if (!group_col %in% names(df)) stop("Grouping column '", group_col, "' not found in res_anno")
  if (!substance_col %in% names(df)) stop("Substance column '", substance_col, "' not found in res_anno")

  # summarize per pathway and per substance (level = substance)
  tbl <- df %>%
    dplyr::group_by_at(c(group_col, substance_col)) %>%
    dplyr::summarize(
      nb_of_items = dplyr::n(),
      firstquartile = as.numeric(stats::quantile(BMD.zSD, 0.25, na.rm = TRUE)),
      secondquartile = as.numeric(stats::median(BMD.zSD, na.rm = TRUE)),
      thirdquartile = as.numeric(stats::quantile(BMD.zSD, 0.75, na.rm = TRUE)),
      min_BMD = as.numeric(min(BMD.zSD, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::ungroup()

  # rename grouping column to 'groupby'
  names(tbl)[names(tbl) == group_col] <- "groupby"

  # optionally filter groups with few items
  tbl <- tbl %>% dplyr::filter(nb_of_items >= min_items)

  tbl
}