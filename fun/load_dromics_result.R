#' Load DRomics Result
#'
#' This function attempts to load a DRomics result file for a given substance. 
#' It searches for the file in multiple predefined paths and loads the first 
#' file it finds. If no file is found, an error is raised.
#'
#' @param substance A character string representing the name of the substance 
#' for which the DRomics result file is to be loaded.
#'
#' @return The content of the DRomics result file as an R object, loaded using 
#' `readRDS()`.
#'
#' @details The function searches for the result file in the following paths:
#' \itemize{
#'   \item `<repo_root>/<substance>/results/<substance>_DRomics_results.rds`
#'   \item `<repo_root>/data/<substance>/results/<substance>_DRomics_results.rds`
#'   \item `<repo_root>/<substance>/results/<substance>_DRomics_result.rds`
#'   \item `<repo_root>/data/<substance>/results/<substance>_DRomics_result.rds`
#' }
#' The first file found in the above paths is loaded and returned. If no file 
#' is found, the function stops with an error message listing the paths it 
#' searched.
#'
#' @examples
#' \dontrun{
#' # Load DRomics result for a substance named "example_substance"
#' result <- load_dromics_result("example_substance")
#' }
#'
#' @export

load_dromics_result <- function(substance) {
  candidate_paths <- c(
    file.path(repo_root, substance, "results", paste0(substance, "_DRomics_results.rds")),
    file.path(repo_root, "data", substance, "results", paste0(substance, "_DRomics_results.rds")),
    file.path(repo_root, substance, "results", paste0(substance, "_DRomics_result.rds")),
    file.path(repo_root, "data", substance, "results", paste0(substance, "_DRomics_result.rds"))
  )
  found <- candidate_paths[file.exists(candidate_paths)]
  if (length(found) == 0) stop("Could not find DRomics result for '", substance, "'. Looked for:\n", paste(candidate_paths, collapse = "\n"))
  readRDS(found[1])
}
