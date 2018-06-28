#' use scImpute to impute dropout values in scRNA-seq data
#'
#' @param count_path A character specifying the full path of the raw count matrix;
#' @param infile A character specifying the type of file storing the raw count matrix;
#' can be "csv", "txt", or "rds". The input file shoule have rows representing genes and
#' columns representing cells, with its first row as cell names 
#' and first column as gene names.
#' @param outfile A character specifying the type of file storing the imputed count matrix;
#' can be "csv", "txt", or "rds".
#' @param out_dir A character specifying the full path of the output directory, 
#' which is used to store all intermdediate and final outputs.
#' @param type A character specifying the type of values in the expression matrix. 
#' Can be "count" (default) or "TPM".
#' @param labeled A logical value indicating whether cell type information is available.
#' \code{labels} must be specified if \code{labeled = TRUE}.
#' @param genelen An integer vector giving the length of each gene. 
#' Order must match the gene orders in the expression matrix.
#' \code{genelen} must be specified if \code{type = "count"}.
#' @param drop_thre A number between 0 and 1, 
#' specifying the threshold to determine dropout values.
#' @param Kcluster An integer specifying the number of cell subpopulations. 
#' This parameter can be determined based on prior knowledge or clustering of raw data.
#' \code{Kcluster} is used to determine the candidate neighbors of each cell.
#' @param labels A character vector specifying the cell type of 
#' each column in the raw count matrix. Only needed when \code{labeled = TRUE}.
#' Each cell type should have at least two cells for imputation.
#' @param ncores A integer specifying the number of cores used for parallel computation.
#' @return scImpute returns a vector giving the column indices of outlier cells.
#' It saves the imputed count matrix to scimpute_count.csv, scimpute_count.txt, or scimpute_count.rds 
#' (depending on \code{outfile}) to \code{out_dir}.
#' @export
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom stats complete.cases dgamma dnorm prcomp quantile rgamma rnorm sd uniroot
#' @importFrom kernlab specc
#' @import penalized 
#' @importFrom utils read.csv read.table write.csv write.table
#' @importFrom rsvd rpca
#' @author Wei Vivian Li, \email{liw@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
#' @references Li, W. V., & Li, J. J. (2018). An accurate and robust imputation method 
#' scImpute for single-cell RNA-seq data. \emph{Nature Communications}, 9(1), 997.
scimpute <-
function (count_path, infile = "csv", outfile = "csv", type = "count", out_dir, labeled = FALSE, 
          drop_thre = 0.5, Kcluster = NULL, labels = NULL, genelen = NULL, ncores = 5) 
{   
    if(labeled == TRUE & is.null(labels)){
      stop("'labels' must be specified when 'labeled = TRUE'!")
    }
    if(labeled == FALSE & is.null(Kcluster)){
      stop("'Kcluster' must be specified when 'labeled = FALSE'!")
    }
    if(!(type %in% c("count", "TPM"))){ stop("expression values can be either 'count' or 'TPM'!") }
    if(type == "TPM" & is.null(genelen)){ stop("'genelen' must be specified when type = 'TPM'!") }
  
    # print(drop_thre)
    print("reading in raw count matrix ...")
    dir.create(out_dir, recursive = TRUE)
    count_lnorm = read_count(filetype = infile, path = count_path, out_dir = out_dir, 
                             type = type, genelen = genelen)
    print("reading finished!")
    
    if(labeled == TRUE){
      if(length(labels) != ncol(count_lnorm)){
        stop("number of cells does not match number of labels !")
      }
    }
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    
    print("imputation starts ...")
    if (labeled == FALSE){
      res_imp = imputation_model8(count = count_lnorm, labeled = FALSE, 
                                  point = log10(1.01), drop_thre = drop_thre, 
                                  Kcluster = Kcluster, 
                                  out_dir = out_dir, ncores = ncores)
    }else{
      res_imp = imputation_wlabel_model8(count = count_lnorm, labeled = TRUE, 
                                         cell_labels = labels, point = log10(1.01), 
                                         drop_thre = drop_thre, 
                                         Kcluster = NULL, out_dir = out_dir, 
                                         ncores = ncores)
    }
    count_imp = res_imp$count_imp
    outliers = res_imp$outlier
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    write_count(count_imp, filetype = outfile, out_dir, type = type, genelen = genelen)
    return(outliers)
}



