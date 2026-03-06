# genx-pfoa-zebrafish-tpod-r
# Sensitive Transcriptomic Points of Departure for GenX and PFOA – Analysis Code

This repository provides R code to reproduce the main analyses for the publication  
“Sensitive Transcriptomic Points of Departure for GenX and PFOA: Implications for PFAS Risk Assessment Using Zebrafish Embryos.”  
It covers transcriptomic point of departure (tPOD) derivation, and pathway‑level sensitivity assessment based on zebrafish embryo (OECD TG 236) transcriptomics.

## Structure

- `main.R`: End‑to‑end analysis script.  
  - imports preprocessed count matrix and sample annotation,  
  - performs identification of regulated genes,  
  - derives tPODs and most sensitive pathways,  
  - produces summary tables and plots for PFOA and GenX.
- `fun`: folder with utility functions for package loading, safe reading, data wrangling, and plotting.


## Execution

1. Place the processed count matrix and sample annotation file(s) in the repository (default: `data/` or project root, as expected in `main.R`).  
2. Open an R session in the project root or set the working directory accordingly.  
3. Install required packages once 
4. Run the analysis:
   - interactively by sourcing `main.R`, or  
   - from the command line via  
     `Rscript main.R`

All key analysis steps executed by the paper (tPOD, pathway sensitivity) are triggered from `main.R`.

## Outputs

The script generates, for each substance and analysis step:

- Figures (e.g. volcano plots, tPOD distributions, pathway sensitivity plots) written to `results/` (and subfolders if defined in `main.R`).  
- Tabular summaries (e.g. differential expression results, tPOD tables, pathway‑level metrics) in CSV format in `results/`.

Exact file names and subfolders are documented in comments at the top of `main.R`.

## Data availability

This repository does **not** contain raw sequencing data.  
Input data for the analyses (count matrix and sample metadata) can be obtained from:

- [insert Zenodo / GEO / institutional repository DOI or URL here]  

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

> Wähle hier eine Lizenz (z.B. `MIT` oder `GPL-3.0`) und passe die Zeile entsprechend an.
