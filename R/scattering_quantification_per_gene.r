#' An expression scattering measurement function
#'
#' An expression scattering measurement function computes the level of scattering 
#' for individual genes along the cell ordering
#'
#' @param CDS a Monocle's CellDataSet object
#' @return integer
#' @author MaiChan Lau
#' @export

scattering_quantification_per_gene <- function(CDS = NULL){

  x <- CDS@assayData$exprs
  if (CDS@expressionFamily@vfamily %in% c("zanegbinomialff", "negbinomial",
                                      "poissonff", "quasipoissonff")) {
    expression <- round(x)
    integer_expression <- TRUE
  }
  else if (CDS@expressionFamily@vfamily %in% c("gaussianff")) {
    expression <- x
  }
  else {
    expression <- log10(x)
    integer_expression <- FALSE
  }

  log_exprs <- log10(x)
  log_exprs[!is.finite(log_exprs)] <- log10(CDS@lowerDetectionLimit)

  total.dispersion <- apply(log_exprs, 1, variability_per_gene)
  return(total.dispersion)
}

#' A utility function for scattering_quantification_per_gene
#'
#' A utility function for scattering_quantification_per_gene which computes the degree 
#' of scattering for single gene, whereby the value is computed by summing over the
#' local values of smaller local windows
#'
#' @param logExp a log-scale expression vector of a gene
#' @param min_expr a minimum expression value
#' @param window_size_perct a window size (in % of total no. of cells) used for computing
#' dispersion level
#' @param nonZeroExpr_perct a minimum amount of cells (in % of window size) with non-zero
#' expression, otherwise the associated window will be assigned to 0 disperson value
#' @return integer
#' @author MaiChan Lau
#' @export

variability_per_gene <- function(logExp =  NULL, min_expr = 0.1,
                                 window_size_perct = 0.1, nonZeroExpr_perct = 0.1){

  num_cell<-length(logExp)
  window_size<-ceiling(num_cell*window_size_perct)
  ind.cell.var <- 0
  total.dispersion<-0

  for(i in 1:num_cell){

    data_window<-logExp
    startID <-(i-1)*window_size + 1
    endID <- startID + (window_size)
    if(endID > num_cell) endID <- num_cell

    data_window<-data_window[startID:endID]
    if(sum(data_window > log10(min_expr)) <= nonZeroExpr_perct*window_size) ind.cell.var <- 0

    else {
      window_nonZero <- data_window[data_window > log10(min_expr)]
      ind.cell.var <- diff(range(window_nonZero))/diff(range(logExp))
    }
    total.dispersion<-total.dispersion+ind.cell.var
    if(endID==num_cell) break;
  }

  return(total.dispersion)
}



