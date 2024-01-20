#' Perform multiple test correction on rows or columns of a p-value matrix
#' 
#' @description
#' User chooses if multiple test correction should be performed row-wise or column-wise.
#' Depending on this choice, rows or columns are treated as the p-values from
#' independent hypothesis tests and are corrected using the [stats::p.adjust()] function.
#' 
#' @param pdat matrix of pvalues
#' @param dim dimension along which adjustment should be performed
#' \itemize{
#' \item if \code{dim==1}, each row is treated as a collection of hypothesis tests for which
#'  a set of q-values should be computed
#'  \item if \code{dim==2}, each column is treated as a collection of hypothesis tests for which
#'  a set of q-values should be computed
#' }
#' @param pmethod adjustment method, passed to [stats::padjust(method=pmethod)]
#' 
#' @return matrix of q-values
#' @author Turner Silverthorne
#' @export
matrix_1d_padjust=function(pdat,dim,pmethod){
  qdat=NaN*pdat
  
  #TODO: could write more concisely using apply() but this gave cryptic rann errors
  if (dim==1){ 
    for (ii in c(1:dim(pdat)[1])){
      qdat[ii,] = p.adjust(pdat[ii,],method=pmethod) 
    }
  }
  
  if (dim==2){
    for (ii in c(1:dim(pdat)[2])){
      qdat[,ii] = p.adjust(pdat[,ii],method=pmethod) 
    }
  }
    
  
  return(qdat)
}