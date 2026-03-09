#' @title Retrieve KEGG Pathway Annotations
#' 
#' @description This function queries the KEGG database to retrieve pathway annotations for a specified species. 
#' It maps KEGG pathway IDs to gene IDs and further maps Entrez gene IDs to ENSEMBL gene IDs using the appropriate 
#' annotation database. The resulting annotation can be saved to a file if desired.
#' 
#' @param species A character string specifying the species code for KEGG queries. Default is "dre" (Danio rerio, zebrafish).
#' @param save_path A character string specifying the file path to save the resulting annotation data frame. 
#' If NULL, the data frame is not saved. Default is NULL.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' 
#' @return A data frame containing the KEGG pathway annotations, including pathway IDs, pathway names, Entrez gene IDs, 
#' and ENSEMBL gene IDs.
#' 
#' @details 
#' The function performs the following steps:
#' \itemize{
#'   \item Queries KEGG for pathway and gene-to-pathway mappings for the specified species.
#'   \item Maps Entrez gene IDs to ENSEMBL gene IDs using the appropriate annotation database (e.g., org.Dr.eg.db for zebrafish).
#'   \item Optionally saves the resulting annotation data frame to a specified file path.
#' }
#' 
#' @note The function requires the \code{KEGGREST}, \code{AnnotationDbi}, and species-specific annotation packages 
#' (e.g., \code{org.Dr.eg.db} for zebrafish) to be installed.
#' 
#' @importFrom KEGGREST keggList keggLink
#' @importFrom AnnotationDbi select
#' @importFrom dplyr mutate left_join distinct
#' @importFrom utils write.table
#' 
#' @examples
#' \dontrun{
#' # Retrieve KEGG annotations for zebrafish and save to a file
#' annotations <- get_kegg_annotation(species = "dre", save_path = "zebrafish_kegg_annotations.txt")
#' 
#' # Retrieve KEGG annotations for zebrafish without saving
#' annotations <- get_kegg_annotation(species = "dre")
#' }
#' 
#' @export
get_kegg_annotation <- function(species = "dre", save_path = NULL, verbose = TRUE) {
  # Build KEGG pathway <-> gene mapping and map Entrez IDs to ENSEMBL IDs
  if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("KEGGREST required for KEGG annotation")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("AnnotationDbi required for ID mapping")
  if (!requireNamespace("org.Dr.eg.db", quietly = TRUE)) stop("org.Dr.eg.db required for zebrafish ID mapping")

  if (verbose) message("Querying KEGG for pathways (species = ", species, ")...")
  zebrafish_pathways <- KEGGREST::keggList("pathway", species)
  zebrafish_genes_to_pathways <- KEGGREST::keggLink("pathway", species)

  pathways_df <- data.frame(
    pathway_id = names(zebrafish_pathways),
    pathway_name = unname(zebrafish_pathways),
    stringsAsFactors = FALSE
  )

  genes_to_pathways_df <- data.frame(
    gene_id = sub(paste0(species, ":"), "", names(zebrafish_genes_to_pathways)),
    pathway_id = sub("path:", "", zebrafish_genes_to_pathways),
    stringsAsFactors = FALSE
  )

  annotation_df <- merge(genes_to_pathways_df, pathways_df, by = "pathway_id", all.x = TRUE)
  annotation_df$pathway_name <- gsub(" - Danio rerio \\(zebrafish\\)", "", annotation_df$pathway_name)

  # Map Entrez -> ENSEMBL using org.Dr.eg.db
  if (verbose) message("Mapping Entrez IDs to ENSEMBL IDs via org.Dr.eg.db...")
  id_mapping <- AnnotationDbi::select(org.Dr.eg.db,
                                      keys = as.character(annotation_df$gene_id),
                                      columns = c("ENSEMBL", "ENTREZID"),
                                      keytype = "ENTREZID")

  annotation_df <- annotation_df %>% dplyr::mutate(gene_id = as.character(gene_id))
  annot_with_ensembl <- annotation_df %>%
    dplyr::left_join(id_mapping, by = c("gene_id" = "ENTREZID")) %>%
    dplyr::distinct(gene_id, ENSEMBL, .keep_all = TRUE)

  if (!is.null(save_path)) {
    utils::write.table(annot_with_ensembl, file = save_path, sep = "\t", row.names = FALSE, quote = FALSE)
    if (verbose) message("Saved KEGG annotation to ", save_path)
  }

  annot_with_ensembl
}