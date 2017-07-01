#' use SCimpute to impute dropout values in scRNA-seq data
#'
#' @param count_path A character specifying the full path of the raw count matrix;
#' @param infile A character specifying the type of file storing the raw count matrix;
#' can be either "csv" or "txt". The input file shoule have rows representing genes and
#' columns representing cells, with its first row as cell names 
#' and first column as gene names.
#' @param outfile A character specifying the type of file storing the imputed count matrix;
#' can be either "csv" or "txt".
#' @param out_dir A character specifying the full path of the output directory, 
#' which is used to store all intermdediate and final outputs.
#' @param drop_thre A number between 0 and 1, 
#' specifying the threshold to determine dropout values.
#' @param celltype A logical value indicating whether cell type information is available.
#' \code{labels} must be specified if \code{celltype = TRUE}.
#' @param labels A character vector specifying the cell type of 
#' each column in the raw count matrix. Only needed when \code{celltype = TRUE}.
#' Each cell type should have at least two cells for imputation.
#' @param ncores A integer specifying the number of cores used for parallel computation.
#' @return Save the imputed count matrix to SCimpute.csv or SCimpute.txt 
#' (depending on \code{outfile}) to \code{out_dir}.
#' @export
#' @import parallel
#' @import glmnet
#' @import stats
#' @import utils
scimpute <-
function (count_path, infile = "csv", outfile = "csv", out_dir, 
    drop_thre = 0.5, celltype = FALSE, labels = NULL, ncores = 5) 
{   
    if(celltype == TRUE & is.null(labels)){
      print("'labels' must be specified when 'celltype = TRUE'!"); stop()
    }
    # print(drop_thre)
    print("reading in raw count matrix ...")
    count_lnorm = read_count(filetype = infile, path = count_path, out_dir = out_dir)
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    print("estimating mixture models ...")
    get_mix_parameters(count = count_lnorm, point = log10(1.01), 
        path = paste0(out_dir, "parslist.rds"), ncores = ncores)
    parslist = readRDS(paste0(out_dir, "parslist.rds"))
    
    print("imputing dropout values ...")
    if (celltype == FALSE){
      count_imp = imputation_model1(count = count_lnorm, point = log10(1.01), 
                                    parslist, drop_thre = drop_thre, method = 2, ncores = ncores)
    }else{
      count_imp = imputation_model1_bytype(count = count_lnorm, labels, point = log10(1.01), parslist, 
                                           drop_thre = drop_thre, method = 2, ncores = ncores)
    }
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    write_count(count_imp, filetype = outfile, out_dir)
    return(0)
}



#' quick re-run of SCimpute with a different \code{drop_thre}
#'
#' @param count_path A character specifying the full path of the raw count matrix;
#' @param infile A character specifying the type of file storing the raw count matrix;
#' can be either "csv" or "txt". The input file shoule have rows representing genes and
#' columns representing cells, with its first row as cell names 
#' and first column as gene names.
#' @param outfile A character specifying the type of file storing the imputed count matrix;
#' can be either "csv" or "txt".
#' @param out_dir A character specifying the full path of the output directory, 
#' which is used to store all intermdediate and final outputs.
#' @param drop_thre A number between 0 and 1, 
#' specifying the threshold to determine dropout values.
#' @param celltype A logical value indicating whether cell type information is available.
#' \code{labels} must be specified if \code{celltype = TRUE}.
#' @param labels A character vector specifying the cell type of 
#' each column in the raw count matrix. Only needed when \code{celltype = TRUE}.
#' Each cell type should have at least two cells for imputation.
#' @param ncores A integer specifying the number of cores used for parallel computation.
#' @return Save the imputed count matrix to SCimpute.csv or SCimpute.txt 
#' (depending on \code{outfile}) to \code{out_dir}.
#' @export
#' @import parallel
#' @import glmnet
#' @import stats
#' @import utils
scimpute_quick <-
  function (count_path, infile = "csv", outfile = "csv", out_dir, 
            drop_thre = 0.5, celltype = FALSE, labels = NULL, ncores = 5) 
  {
    if(celltype == TRUE & is.null(labels)){
      print("'labels' must be specified when 'celltype = TRUE'!"); stop()
    }
    # print(drop_thre)
    print("reading in raw count matrix ...")
    count_lnorm = read_count(filetype = infile, path = count_path, out_dir = out_dir)
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    parslist = readRDS(paste0(out_dir, "parslist.rds"))
    print("imputing dropout values ...")
    if (celltype == FALSE){
      count_imp = imputation_model1(count = count_lnorm, point = log10(1.01), 
                                    parslist, drop_thre = drop_thre, method = 2, ncores = ncores)
    }else{
      count_imp = imputation_model1_bytype(count = count_lnorm, labels = labels, point = log10(1.01), parslist, 
                                           drop_thre = drop_thre, method = 2, ncores = ncores)
    }
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    write_count(count_imp, filetype = outfile, out_dir)
    return(0)
  }
