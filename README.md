scImpute: a statistical method for accurate and robust imputation of scRNA-seq data
================
Wei Vivian Li, Jingyi Jessica Li
2017-05-24

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

The emerging single cell RNA sequencing (scRNA-seq) technologies enable the investigation of transcriptomic landscape at single-cell resolution. However, scRNA-seq analysis is complicated by the excess of zero or near zero counts in the data, which are the so-called dropouts due to low amounts of mRNA within each individual cell. Consequently, downstream analysis of scRNA-seq woule be severely biased if the dropout events are not properly corrected. `scImpute` is developed to accurately and robustly impute the dropout values in scRNA-seq data.

`scImpute` can be applied to raw read count matrix before the users perform downstream analyses such as

-   dimension reduction of scRNA-seq data
-   normalization of scRNA-seq data
-   clustering of cell populations
-   differential gene expression analysis
-   time-series analysis of gene expression dynamics

Any suggestions on the package are welcome!

Installation
------------

The package is not on CRAN yet. For installation please use the following codes in `R`

``` r
install.packages("devtools")
library(devtools)

install_github("Vivianstats/scImpute")
```

Quick start
-----------

`scImpute` can be easily incorporated into existing pipeline of scRNA-seq analysis. Its only input is the raw count matrix with rows representing genes and columns representing cells. It will output an imputed count matrix with the same dimension. In the simplest case, the imputation task can be done with one single function `scimpute`:

``` r
scimpute(# full path to raw count matrix
         count_path = system.file("extdata", "raw_count.csv", package = "scImpute"), 
         infile = "csv",           # format of input file
         outfile = "csv",          # format of output file
         out_dir = "./",           # full path to output directory
         drop_thre = 0.5,          # threshold set on dropout probability
         ncores = 10)              # number of cores used in parallel computation
```

This function will create a new file `scimpute_count.csv` in `out_dir` to store the imputed count matrix. For detailed usage, please refer to the package manual or vignette.
