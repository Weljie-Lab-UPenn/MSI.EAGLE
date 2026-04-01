<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSI.EAGLE

MSI.EAGLE is an interactive platform for mass spectrometry imaging (MSI)
analysis, designed for both routine MSI workflows and advanced multimodal
workflows. It combines MSI preprocessing, visualization, embedding, and
statistics with optional polygon-informed single-cell analysis and
registration modules in one reproducible interface.

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

### Prerequisites

Before installing MSI.EAGLE, you need to install R and R Studio: https://rstudio-education.github.io/hopr/starting.html

Then the required Bioconductor packages and Git.

### Step 1: Install Required R Packages

```r
# Install remotes package
install.packages("remotes")

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c("BiocParallel", "Cardinal"))

# Optional: Install CardinalWorkflows for demo data
BiocManager::install("CardinalWorkflows")
```

### Step 2: Install Git

If you don't have Git installed, download and install it from: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

### Step 3: Install MSI.EAGLE

```r
#Main branch (stable):
remotes::install_git("https://github.com/Weljie-Lab-UPenn/MSI.EAGLE")
```

OR

```r
#Development branch (recent updates)
remotes::install_git("https://github.com/Weljie-Lab-UPenn/MSI.EAGLE", ref = 'dev')
```

## Usage

### Basic Setup

In R Studio:
1) Create a new R script, and copy and paste the start code below into it.
2) Add the file path to the folder containing your .imzML and .ibd raw data files (rawd), and the path to the folder where your results will be saved (wd).
3) Then, click "Source with Echo" to run the application.

```r
#MSI.EAGLE Start

library(MSI.EAGLE)

# Set path to raw data files
rawd = 'path/to/rawfilesdirectory'

# Set working directory (can be same as raw data directory)
wd = 'path/to/workingdirectory'

# Configure number of CPU cores to use
# Option 1: Use half of available processors
ncores <- as.integer(parallel::detectCores()/2)

# Option 2: Use all but 2 processors (recommended)
ncores <- as.integer(parallel::detectCores()) - 2

# Launch the application
MSI.EAGLE()
```
## Documentation
For detailed instructions, tutorials, and examples, see our Manuals and Vignettes.

## Requirements

- R version 3.6 or higher
- Bioconductor packages: Cardinal, BiocParallel
- Git (for installation)
- Sufficient RAM for mass spectrometry imaging data processing

## Support

For issues and questions, please use the GitHub Issues page for this repository.
