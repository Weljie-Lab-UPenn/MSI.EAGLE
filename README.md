
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSI.EAGLE

<!-- badges: start -->
<!-- badges: end -->

The goal of MSI.EAGLE is to create an interactive user app for
processing of mass spec imaging data. It uses many of the processing
features from Cardinal, and adds a number of analysis options.

## Installation

You can install the development version of MSI.EAGLE as follows:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
Install using: 

Step 1*: Install Bioconductor packages

install.packages("remotes")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("BiocParallel", "Cardinal"))
```

    Step 2: Install Git

<https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>


    Step 3: Install MSI.EAGLE

    remotes::install_git("http://weljie.myds.me:30003/aalim/DESI_Shiny_Processing_script.git", branch="modules")


    ## Example

    This is a basic example which shows you how to solve a common problem:


    ``` r
    ## basic example code
    # wd = ('path/to/workingdirectory)
    # rawd = ('path/to/rawfilesdirectory') # rawd = wd
    # ncores = as.integer(parallel::detectCores()/2) # use 1/2 of the available processors
    # ncores = as.integer(parallel::detectCores())-2 # us all but 2 available processors
    # library(MSI.EAGLE)
    # MSI.EAGLE()

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
