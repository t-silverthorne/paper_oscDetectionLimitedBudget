#' Helper function computes reduced Fisher information matrix for cosinor model
#' 
#' Here, reduced FIM correpsonds to de-meaned model. 
#' @export
getReducedFIM <- function(t,param){
  freq   = param[['freq']];
  
  Y = matrix(c(cos(2*pi*freq*t),
           sin(2*pi*freq*t)),
           nrow=2,byrow=T)
  M=Y %*% t(Y)
  if (any(is.na(M))){
    stop('NA value in FIM')
  }else{
    return(M)
  }
}