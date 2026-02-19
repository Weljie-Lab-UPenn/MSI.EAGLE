<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSI.EAGLE

<!-- badges: start -->
<!-- badges: end -->

MSI.EAGLE is an interactive platform for mass spectrometry imaging (MSI) analysis, designed for both routine MSI workflows and advanced multimodal workflows. It combines MSI preprocessing, visualization, embedding, and statistics with optional polygon-informed single-cell analysis and registration modules in one reproducible interface.

## Why MSI.EAGLE

- General MSI workflow by default: load imzML/RDS data, preprocess, visualize ions, and run clustering/statistics without requiring histology or polygons.
- Optional single-cell workflows: import cell polygons from microscopy, map them to MSI pixels, and store cell-level fields in `pData`.
- Spatially aware UMAP and clustering: explore embedding space and tissue space together.
- Multimodal registration: align histology images, cluster maps, and polygons to MSI with tunable, saveable transforms.
- Quantitative registration diagnostics: score candidate transforms and review fit statistics before applying.
- Reproducible export: save/load registration parameters and write updated MSI datasets with mapped metadata.

## Core Workflows

1. Load imzML or RDS MSI datasets.
2. Preprocess and peak-pick using Cardinal-based methods.
3. Generate UMAP/clustering and derive phenotypes.
4. Optionally register histology, cluster, and polygon overlays and map polygon/cell annotations into `pData`.
5. Run group statistics at pixel or sample level and export results.

## Installation

``` r
install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("BiocParallel", "Cardinal"))

# Optional demo data:
BiocManager::install("CardinalWorkflows")

# MSI.EAGLE:
remotes::install_git("https://github.com/Weljie-Lab-UPenn/MSI.EAGLE")
```

## Quick Start

``` r
# optional globals used by the app
rawd <- "path/to/rawfilesdirectory"
wd <- "path/to/workingdirectory" # can be same as rawd

# use all but 2 available processors (minimum 1)
ncores <- max(1L, as.integer(parallel::detectCores()) - 2L)

library(MSI.EAGLE)
# MSI.EAGLE() # uncomment to start the app
```

## Access to Private Repositories

If needed, store your GitHub Personal Access Token (PAT) in `.Renviron`:

``` r
usethis::edit_r_environ()
```

Add:

``` text
GITHUB_PAT=your_personal_access_token
```

Then restart R.
