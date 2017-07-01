scImpute: accurate and robust imputation of scRNA-seq data
================
Wei Vivian Li, Jingyi Jessica Li
2017-07-01

<!-- README.md is generated from README.Rmd. Please edit that file -->
News
----

> 2017/07/01:

-   Version 0.0.2 is released!
-   This version speeds up the first step in `scImpute` and program now completes in a few seconds when applied to a dataset with 10,000 genes and 100 cells (using single core).

Introduction
------------

`scImpute` is developed to accurately and robustly impute the dropout values in scRNA-seq data. `scImpute` can be applied to raw read count matrix before the users perform downstream analyses such as

-   dimension reduction of scRNA-seq data
-   normalization of scRNA-seq data
-   clustering of cell populations
-   differential gene expression analysis
-   time-series analysis of gene expression dynamics

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/scImpute/issues). For suggestions and comments on the method, please contact Wei (<liw@ucla.edu>) or Dr. Jessica Li (<jli@stat.ucla.edu>).

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

This function will create a new file `scimpute_count.csv` in `out_dir` to store the imputed count matrix. For detailed usage, please refer to the package [manual](https://github.com/Vivianstats/scImpute/blob/master/inst/docs/) or [vignette](https://github.com/Vivianstats/scImpute/blob/master/vignettes/scImpute-vignette.Rmd).
