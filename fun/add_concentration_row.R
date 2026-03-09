#' Add concentration row for DRomics inputs
#'
#' Prepends the concentration vector to the count matrix and returns a DRomics-
#' ready data frame with an `ensembl_gene_id` column. Handles two common
#' coldata layouts (new multi-column layout and legacy `mConc`/`Conc` columns).
#'
#' @param col data.frame Sample metadata (coldata).
#' @param matrix data.frame Count matrix (rows = genes, cols = samples).
#' @param conc_type character Concentration column type(s) in the new layout.
#' @param conc_col character Optional explicit concentration column name.
#'   If provided, this takes precedence over auto-detection in the legacy layout.
#'
#' @return A data.frame (or list of data.frames for the new layout) containing
#'   `ensembl_gene_id` followed by sample columns with concentrations in the
#'   first row.
#' @export
add_concentration_row <- function(col = coldata,
                                  matrix = countmatrix,
                                  conc_type = c("nominal", "initial", "measured"),
                                  conc_col = NULL) {
  # Prüfe auf neue Struktur
  neue_spalten <- c("Conc_nominal", "Conc_initial", "Conc_meanmeasured")
  hat_neue_struktur <- all(neue_spalten %in% colnames(col))
  
  if (hat_neue_struktur) {
    # Mapping von conc_type zu Spaltennamen
    type_map <- list(
      nominal  = "Conc_nominal",
      initial  = "Conc_initial",
      measured = "Conc_meanmeasured"
    )
    
    # Prüfe, ob rownames vorhanden sind, sonst verwende die erste Spalte
    if (!is.null(rownames(col)) && all(rownames(col) != "")) {
      sample_names <- rownames(col)
    } else {
      sample_names <- as.character(col[[1]])
    }
    
    result <- list()
    
    for (typ in conc_type) {
      colname <- type_map[[typ]]
      if (!is.null(colname) && colname %in% colnames(col)) {
        conc <- setNames(col[[colname]], sample_names)
        assign("conc", conc, envir = .GlobalEnv)
        
        new_row <- data.frame(matrix(unlist(conc[colnames(matrix)]), nrow = 1))
        colnames(new_row) <- colnames(matrix)
        
        new_matrix <- rbind(new_row, matrix)
        rownames(new_matrix)[1] <- "ensembl_gene_id"
        
        CM <- suppressWarnings(
          new_matrix %>% rownames_to_column(var = "ensembl_gene_id"))
        
        result[[typ]] <- CM
        write.table(CM, file = paste0(substance,"_", typ, "_DRmatrix.txt"), row.names = F, col.names = F) 
        
      } else {
        warning(paste("Spalte für Typ", typ, "nicht gefunden!"))
      }
    }
    return(result)
  } else {
    # Alte Funktionalität
  mConc_col <- if (!is.null(conc_col)) conc_col else grep("mConc", colnames(col), value = TRUE)
    
    if (!is.null(rownames(col)) && all(rownames(col) != "")) {
      sample_names <- rownames(col)
    } else {
      sample_names <- as.character(col[[1]])
    }
    
    if (length(mConc_col) > 0 && any(!is.na(col[[mConc_col]]))) {
      conc <- setNames(col[[mConc_col]], sample_names)
    } else {
      conc <- setNames(col$Conc, sample_names)
    }
    
    assign("conc", conc, envir = .GlobalEnv)
    
    new_row <- data.frame(matrix(unlist(conc[colnames(matrix)]), nrow = 1))
    colnames(new_row) <- colnames(matrix)
    
    countmatrix2 <- rbind(new_row, matrix)
    rownames(countmatrix2)[1] <- "ensembl_gene_id"
    
    CM <- suppressWarnings(
      countmatrix2 %>% rownames_to_column(var = "ensembl_gene_id"))
    write.table(CM, file = paste0(substance, "_DRmatrix.txt"), row.names = F, col.names = F) 
    
    return(CM)
  }
}
