#'Reduced Fisher information matrix for cosinor model
#' 
#' @param t vector of measurement times
#' @param$freq frequency of oscillation
#' 
#' @return Fisher information matrix corresponding to de-meaned cosinor model. In 
#' the paper, we refer to this as the reduced Fisher information matrix.
#' 
#' @author Turner Silverthorne 
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

#'Helper function returns smallest eigenvalue
#' 
#' @description
#' Wrapper for computing smallest eigenvalue of a matrix
#' 
#' @param Mat a matrix
#' @param is_symmetric boolean, is the matrix symmetric
#' 
#' @return smallest eigenvalue of mat, computed using [base::eigen] 
#' 
#' @author Turner Silverthorne 
#' @export
getMinEig <- function(Mat,is_symmetric=F){
  return(base::eigen(Mat,only.values=T,symmetric = is_symmetric) %>% {.$values} %>%  min())
}

#'Frequency derivative of smallest FIM eigenvalue
#' @description
#' computes derivative of minimum eigenvalue of reduced Fisher information matrix 
#' with respect to the frequency parameter
#' 
#' @param mt vector of measurement times
#' @param freq frequency of signal
#' 
#' @return derivative of minimum eigenvalue of reduced FIM with respect to frequency
#' 
#' @author Turner Silverthorne 
#' @export
getDeigDfreq=function(mt,freq){
  cvec=cos(2*pi*freq*mt)
  svec=sin(2*pi*freq*mt)
  
  M11=t(cvec)%*%cvec
  M22=t(svec)%*%svec
  M12=t(cvec)%*%svec
  
  M  = matrix(c(M11,M12,M12,M22),nrow=2)
  M
  eigs = eigen(M,symmetric = T)
  vmin = eigs$vectors[,which.min(eigs$values)]
  
  dM22     = sum(4*pi*mt*cvec*svec)
  dM12     = sum(2*pi*mt*(2*cvec*cvec-1))
  dMdfreq = matrix(c(-dM22,dM12,dM12,dM22),nrow=2,byrow=T)
  
  return(t(vmin)%*%dMdfreq%*%vmin/(t(vmin)%*%vmin))
}
