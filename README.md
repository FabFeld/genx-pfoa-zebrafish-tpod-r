# **Sensitive Transcriptomic Points of Departure for GenX and PFOA – Analysis Code**

This repository provides R code to reproduce the main analyses for the publication\
“Sensitive Transcriptomic Points of Departure for GenX and PFOA: Implications for PFAS Risk Assessment Using Zebrafish Embryos.”\
It covers transcriptomic point of departure (tPOD) derivation, and pathway‑level sensitivity assessment based on zebrafish embryo (OECD TG 236) transcriptomics.

## Repository layout

```         
publication_repro/            # orchestration scripts and helper functions
    fun/                       # helper functions (plotting, annotation, summaries)
    main.R                     # per-substance DRomics driver (runs multiple substances)
    substance_comparison.R     # stand-alone comparison runner
    comparison_plots/          # outputs from the comparison stage
GenX/                        # example substance folder (replace with your own data)
    GenX_coldata.csv
    GenX_CountMatrix.txt
    results/                   # created by main.R
PFOA/                        # second example substance
    PFOA_coldata.csv
    PFOA_CountMatrix.txt
    results/
results/                     # legacy comparison outputs (kept for compatibility)
    substance_comparison/
```

## Overview

`publication_repro/main.R` is the single entry point for the full pipeline. It loops over the configured substance folders, loads metadata/count matrices, annotates Ensembl IDs, runs DRomics selection → fitting → BMD estimation, and saves outputs in `<Substance>/results/`. The same script then generates comparison plots by sourcing the documented helpers in `publication_repro/fun/`.

## Dependencies

Missing packages are installed automatically via `setup_packages()`. The pipeline relies on:

-   **Annotation**: `biomaRt`, `AnnotationDbi`, `org.Dr.eg.db`, `KEGGREST`
-   **Data manipulation**: `dplyr`, `tibble`, `tidyr`, `purrr`, `glue`, `reshape2`, `scales`
-   **Visualization**: `ggplot2`, `ggsci`, `grid`, `gridExtra`, `gridtext`
-   **Modeling/performance**: `DRomics`, `parallel`

## Per-substance workflow

1.  Place `<Substance>_coldata.csv` and `<Substance>_CountMatrix.txt|.csv` inside `<repo_root>/<Substance>/`.
2.  Ensure the metadata contains a concentration column named `mConc_ug.L` (or update `main.R::conc_col`).
3.  Run from the repository root:

``` powershell
Rscript publication_repro/main.R
```

Outputs per substance:

-   `<Substance>_DRomics_results.rds`: saved DRomics object.
-   `<Substance>_summary_table.csv`: combined summary with `method`, `BMD.zSD`, and annotation.
-   `<Substance>_drc_<method>.png`: dose–response plots from `universal_drcplot()`.
-   `<Substance>_session_info.txt`: session info for reproducibility.

## Cross-substance comparison

After the per-substance loop, `main.R` performs the comparison stage:

-   Loads results with `load_dromics_result()`.
-   Builds BMD density plots with `generate_density_plots()`.
-   Retrieves KEGG annotations using `get_kegg_annotation()`.
-   Generates pathway sensitivity summaries and plots with `generate_sensitivity_plot()`.

Comparison outputs are stored in `publication_repro/comparison_plots/`:

-   `density_comparison.pdf`
-   `sensitivity_plot_table.csv`
-   `sensitivity_plots.pdf`
-   `session_info.txt`

## Helper catalog (`publication_repro/fun/`)

-   `add_concentration_row()` – prepares DRomics input by adding the concentration row and writing `*_DRmatrix.txt`.
-   `load_dromics_result()` – robust loader for `<Substance>_DRomics_results.rds`.
-   `generate_density_plots()` / `dens_plot*()` – publication-ready density plots.
-   `build_sensitivity_table()` / `generate_sensitivity_plot()` – pathway-level summaries and sensitivity plots.
-   `get_kegg_annotation()` – caches KEGG mapping (`kegg_annotation.tsv`).
-   `universal_drcplot()` – standard dose–response plot wrapper for DRomics fits.

## Notes

-   `main.R` automatically sources all `.R` files in `publication_repro/fun/`.
-   The script overwrites comparison outputs; archive `comparison_plots/` before re-running if needed.
-   Add new substances by extending `substance_dirs` and providing the correct input files.
-   Functionality was validated only for the documented R version and package versions in `session_info.txt`; results may differ with other versions.

## Data availability

This repository does **not** contain raw sequencing data.  
Input data for the analyses (count matrix and sample metadata) can be obtained from:

- 10.5281/zenodo.18922913

Please download the data and place the files as described in the *Execution* section before running `main.R`.

## Citation

If you use this code, please cite:

- The associated article:  
  *Sensitive Transcriptomic Points of Departure for GenX and PFOA: Implications for PFAS Risk Assessment Using Zebrafish Embryos*  
  [insert journal, year, volume, pages, DOI]

- (Optional) The code repository / Zenodo record if you mint a separate DOI:  
  [insert Zenodo DOI for the code deposition]

## License

This code is distributed under the [MIT/GPL‑3.0] license.  
See the `LICENSE` file for details.

