#!/usr/bin/env Rscript
################################################################################
# publication_repro/main.R
#
# Purpose:
#   Reproducible DRomics pipeline that reads per-substance count matrices and
#   coldata, annotates genes (Ensembl/KEGG), runs DRomics selection/fitting/BMD
#   analyses, and writes results, summary tables and comparison plots.
#
# Usage (interactive):
#   - Open `publication_repro/main.R` in RStudio and run the script (interactive
#     mode will set the working directory to the script folder).
#
# Usage (CLI):
#   Rscript main.R
#
# Inputs (per substance):
#   - <repo_root>/publication_repro/<Substance>/: contains
#       - <Substance>_coldata.csv (sample metadata with concentration column)
#       - <Substance>_CountMatrix.txt (rows = genes, cols = samples)
#
# Outputs:
#   - <repo_root>/publication_repro/<Substance>/results/ : DRomics outputs and plots
#   - <repo_root>/publication_repro/comparison_plots/ : cross-substance comparison CSV/PDF
#
# Dependencies:
#   biomaRt, dplyr, tibble, tidyr, scales, ggplot2, gridExtra, gridtext, grid,
#   reshape2, KEGGREST, org.Dr.eg.db, AnnotationDbi, parallel, DRomics, ggsci
#
# Author: (auto) Positron assistant edits — adapt metadata as needed
# License: CC-BY or your preferred license (not set here)
################################################################################


# ---- Packages -----------------------------------------------------------------
# Minimal fallback for setup_packages if helpers not yet loaded
if (!exists("setup_packages")) {
  setup_packages <- function(pkgs) {
    missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(missing_pkgs) > 0) {
      install.packages(missing_pkgs, dependencies = TRUE, repos = "https://cloud.r-project.org")
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  }
}

# List of required R packages for the script:
# - biomaRt: Provides access to BioMart databases for biological data retrieval (essential).
# - dplyr: A grammar of data manipulation, providing a consistent set of verbs for data wrangling (essential).
# - tibble: A modern reimagining of data frames, making them easier to work with (essential).
# - tidyr: Tools for tidying data, ensuring it is in a consistent format (essential).
# - scales: Provides tools for scaling and formatting data for visualization (non-essential, but useful for plots).
# - ggplot2: A system for declaratively creating graphics based on the Grammar of Graphics (essential for visualization).
# - gridExtra: Provides functions to arrange multiple grid-based plots on a page (non-essential, but useful for complex layouts).
# - gridtext: Enables rich text rendering in grid graphics (non-essential, but useful for enhanced text formatting).
# - grid: Base R package for creating and manipulating graphical objects (essential for low-level graphics).
# - reshape2: Tools for reshaping data between wide and long formats (non-essential, but useful for data transformation).
# - KEGGREST: Provides programmatic access to KEGG databases for pathway and gene data (essential for pathway analysis).
# - org.Dr.eg.db: Annotation package for zebrafish (Danio rerio) gene data (essential for gene annotation).
# - AnnotationDbi: Provides tools for handling and querying annotation data (essential for gene annotation).
# - parallel: Base R package for parallel computing (non-essential, but useful for performance optimization).
# - DRomics: A package for dose-response analysis of omics data (essential for the specific analysis in this script).
# - ggsci: Provides scientific journal and sci-fi themed color palettes for ggplot2 (non-essential, but useful for aesthetics).

# Required packages for the pipeline
required_pkgs <- c(
  "biomaRt", "dplyr", "tibble", "tidyr", "scales", "ggplot2", "gridExtra",
  "gridtext", "grid", "reshape2", "KEGGREST", "org.Dr.eg.db", "AnnotationDbi",
  "parallel", "DRomics", "ggsci"
)
# Configure packages
setup_packages(required_pkgs)

# ---- Paths & working directory ------------------------------------------------
# Resolve script location and set repository root (parent of publication_repro)
repo_root <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_entry <- grep("^--file=", args)
  if (length(file_entry) > 0) {
    normalizePath(sub("^--file=", "", args[file_entry[1]]))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    normalizePath(sys.frames()[[1]]$ofile)
  } else {
    normalizePath(getwd())
  }
}

# If running in RStudio interactive editing session, prefer the active document's folder as working dir
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  doc_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
  if (nzchar(doc_path)) {
    wd_dir <- dirname(doc_path)
    setwd(wd_dir)
  } else {
    setwd(normalizePath(repo_root, winslash = "\\", mustWork = FALSE))
  }
} else {
  # default: use repo_root as working directory
  setwd(normalizePath(repo_root, winslash = "\\", mustWork = FALSE))
}



# ---- Helper sourcing ----------------------------------------------------------
# Load helper utilities from publication_repro/fun
fun_dir <- file.path(repo_root, "fun") |> normalizePath()
if (dir.exists(fun_dir)) {
  fun_files <- list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
  for (fun_file in fun_files) {
    try(source(fun_file), silent = TRUE)
  }
} else {
  warning("Helper functions directory not found: ", fun_dir, "\n",
          "Please ensure you have the publication_repro repository with the 'fun' folder containing helper scripts.")
}



# ---- Analysis configuration ---------------------------------------------------
# Determine substance directories to process
substance_dirs <- c("GenX", "PFOA") # Adjust as needed


# BMD parameters
z <- 2 # number of standard deviations for BMD calculation (zSD)
item_no <- 36 # number of items to plot for dose-response curves
niter <- 1 # number of bootstrap iterations for BMD confidence intervals
select_methods <- c("quadratic") # At least one of c("quadratic", "linear", "ANOVA")
FDR <- 0.05 # False discovery rate threshold for item selection



for (substance in substance_dirs) {
  # ---- Input files ------------------------------------------------------------
  working_dir <- file.path(repo_root, substance) |> normalizePath()
  if (!dir.exists(working_dir)) {
    dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)
    message("Created substance directory: ", working_dir, " — please add your coldata/countmatrix files there.")
  }

  coldata_files <- list.files(
    path = working_dir,
    pattern = ".*[Cc]ol[Dd]ata.*\\.csv$",
    recursive = FALSE,
    full.names = TRUE
  )
  if (length(coldata_files) == 0) {
    stop("No coldata CSV found in ", working_dir)
  }
  coldata <- read.csv2(coldata_files[1], row.names = 1)

  unit <- if ("Unit" %in% colnames(coldata)) unique(coldata$Unit)[1] else expression("µg/L")

  count_files <- list.files(
    path = working_dir,
    pattern = ".*[Cc]ount[Mm]atrix.*\\.(txt|csv)$",
    recursive = FALSE,
    full.names = TRUE
  )
  if (length(count_files) == 0) {
    stop("No count matrix file found in ", working_dir)
  }
  countmatrix <- read.delim(count_files[1], header = TRUE, row.names = 1, check.names = FALSE) |>
    as.data.frame() |>
    dplyr::select(rownames(coldata))

  if (!all(rownames(coldata) %in% colnames(countmatrix))) {
    stop("Sample names in coldata do not match column names in countmatrix.")
  }

  # ---- Annotation -------------------------------------------------------------
  if (!exists("GeneAnno")) {
    annot <- biomaRt::useEnsembl(
      biomart = "ensembl", dataset = "drerio_gene_ensembl",
      host = "https://www.ensembl.org"
    )
    GeneAnno <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description", "entrezgene_id"),
      filters = "ensembl_gene_id",
      values = rownames(countmatrix),
      mart = annot
    ) %>%
      distinct(ensembl_gene_id, .keep_all = TRUE) %>%
      transmute(
        ensembl_gene_id,
        SYMBOL = coalesce(external_gene_name, ensembl_gene_id),
        description,
        entrezgene_id
      )
  }

  # ---- DRomics inputs ---------------------------------------------------------
  conc_col <- "mConc_ug.L"
  if (!conc_col %in% colnames(coldata)) {
    stop("Concentration column missing: ", conc_col)
  }
  conc <- as.numeric(coldata[[conc_col]])
  CM <- add_concentration_row(col = coldata, matrix = countmatrix, conc_col = conc_col)
  CM <- rownames_to_column(CM, var = "ensembl_gene_id")
  o <- RNAseqdata(CM, backgrounddose = min(conc), check = TRUE, transfo.method = "rlog", round.counts = TRUE)

  # ---- DRomics selection, fit, BMD --------------------------------------------
  results <- list()
  for (method in select_methods) {
    message("Processing method: ", method)
    s <- itemselect(o, select.method = method, FDR = FDR)
    if (length(s$selectindex) == 0) next
    f <- drcfit(itemselect = s, progressbar = TRUE)
    f$fitres <- f$fitres %>%
      left_join(GeneAnno, by = c("id" = "ensembl_gene_id")) %>%
      mutate(SYMBOL = coalesce(SYMBOL, id)) %>%
      dplyr::select(id, SYMBOL, everything())
    r <- bmdcalc(f, z = z, x = 10)
    b <- bmdboot(r, niter = niter)
    keep_idx <- f$fitres$irow
    f$omicdata <- list(
      data = f$omicdata$data[keep_idx, , drop = FALSE],
      data.mean = f$omicdata$data.mean[keep_idx, , drop = FALSE],
      raw.counts = f$omicdata$raw.counts[keep_idx, , drop = FALSE],
      item = f$omicdata$item[keep_idx],
      dose = conc
    )
    res <- list(drcfit = f, bmdcalc = r, bmdboot = b, omicdata = f$omicdata)
    res$drcfit$omicdata <- NULL
    res$bmdcalc$omicdata <- NULL
    results[[method]] <- res
  }
  if (length(results) == 0) stop("No DRomics results available.")

  # ---- Outputs ----------------------------------------------------------------
  output_dir <- file.path(working_dir, "results")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()
  for (method in names(results)) {
    res <- results[[method]]
    gg <- universal_drcplot(
      res,
      plot.type = "dose_fitted",
      dose_log_transfo = TRUE,
      gene_no = item_no,
      sorted_by = "maxychange",
      facet_by = "SYMBOL",
      showCI = TRUE
    )
    plots[[method]] <- gg
    ggsave(
      file.path(output_dir, paste0(substance, "_drc_", method, ".png")),
      plot = gg, width = 9, height = 12, dpi = 300
    )
  }

  summary_tables <- lapply(names(results), function(method) {
    res <- results[[method]]
    if (is.null(res)) return(NULL)
    res$drcfit$fitres %>%
      mutate(method = method)
  })

  summary_table <- bind_rows(summary_tables) %>%
    dplyr::select(method, id, SYMBOL, everything())

  write.csv(summary_table, file.path(output_dir, paste0(substance, "_summary_table.csv")), row.names = FALSE)
  saveRDS(results, file.path(output_dir, paste0(substance, "_DRomics_results.rds")))
}






# ---- Comparison stage ---------------------------------------------------------
setwd(repo_root)
# create plots directory for comparison outputs
plots_dir <- file.path(repo_root, "comparison_plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

substances <- c("PFOA", "GenX")

# load results for both substances
results <- setNames(lapply(substances, function(s) { tryCatch(load_dromics_result(s), error = function(e) { warning(e$message); NULL }) }), substances)

if (length(results) == 0) stop("No results found for any substance: check your data/*/results folders.")

generate_density_plots(results, plots_dir)



# ---- Data adaption for functional plots --------------------------------------


# Retrieve KEGG annotation (zebrafish) and keep in `annot_with_ensembl`
if (!exists("annot_with_ensembl")) {
  annot_with_ensembl <- get_kegg_annotation(species = "dre", save_path = file.path(repo_root, "kegg_annotation.tsv"), verbose = TRUE)
}


# bind results from both substances into a combined table for comparison plots (e.g., sensitivity plot by pathway)
res_combined <- bind_rows(
  lapply(names(results), function(substance) {
    res <- results[[substance]]
    method <- "quadratic"
    if (!is.null(res[[method]]) && !is.null(res[[method]]$bmdboot$res)) {
      df <- res[[method]]$bmdboot$res
      df$Substance <- substance
      return(df)
    }
    NULL
  })
)


# Now merge with your results
res_anno <- merge(x = res_combined, y = annot_with_ensembl,
                     by.x = "id", by.y = "ENSEMBL",
                     all.x = TRUE) %>%
  dplyr::filter(!is.na(BMD.zSD) & !is.na(pathway_id) & !is.na(BMD.zSD.lower) & !is.na(BMD.zSD.upper))

# Generate sensitivity table + plot using the documented helper in fun/
sens_out <- generate_sensitivity_plot(
  res_anno,
  output_dir = plots_dir,
  group_col = "pathway_name",
  substance_col = "Substance",
  min_items = 2,
  save_csv = TRUE,
  save_plot = TRUE
)

# expose outputs for downstream use
sensitivity_table <- sens_out$table
sensitivity_plot <- sens_out$plot

message("Substance comparison workflow completed. Outputs written to: ", plots_dir)



# session info for reproducibility
session_info <- sessionInfo()
# Write session info to the comparison plots folder for reproducibility
si_file <- file.path(plots_dir, "session_info.txt")
writeLines(capture.output(session_info), si_file)