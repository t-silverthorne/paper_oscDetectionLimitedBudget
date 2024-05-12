#' Eval worst-case non-centrality parameter
#' @param mt vector of measurement times 
#' @param freqs vector of frequencies at which cost function should be evaluated
#' @param Amp amplitude of signal, only relevant if \code{cfuntype=='power'}
#' @param alpha type I error cutoff, only relevant if \code{cfuntype=='power'}
#' @return minimum value of non-centrality parameter for freq in \code{freqs}
#' @author Turner Silverthorne
#' @export
evalWorstNCP=function(mt,freqs,Amp=1,alpha=.05,method='svd'){
  N=length(mt)  
  cvals = freqs %>% sapply(function(freq){
    X = matrix(c(Amp*cos(2*pi*mt*freq),Amp*sin(2*pi*mt*freq)),ncol=2)
    if (method=='svd'){
      return(svd(X) %>% {.$d} %>% {.^2})
    }else if (method=='eigen'){
      return(eigen(t(X)%*%X) %>% {.$values})
    }else{
      stop('unknown method')
    }
    })
  return(min(cvals))
}
  