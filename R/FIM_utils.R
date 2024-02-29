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

#' Helper function returns smallest eigenvalue
#' @export
getMinEig <- function(Mat,is_symmetric=F){
  return(eigen(Mat,only.values=T,symmetric = is_symmetric) %>% {.$values} %>%  min())
}

#' Derivative of smallest eigenvalue of FIM wrt frequency
deig_dfreq=function(mt,freq){
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
