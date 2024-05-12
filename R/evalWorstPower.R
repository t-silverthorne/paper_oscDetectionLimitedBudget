#' Eval worst-case power 
#' @param mt vector of measurement times 
#' @param freqs vector of frequencies at which cost function should be evaluated
#' @param Amp amplitude of signal, only relevant if \code{cfuntype=='power'}
#' @param alpha type I error cutoff, only relevant if \code{cfuntype=='power'}
#' @return minimum value of non-centrality parameter for freq in \code{freqs}
#' @author Turner Silverthorne
#' @export
evalWorstPower=function(mt,freqs,Amp=1,alpha=.05){
  N   = length(mt)
  f0  = qf(p=1-alpha,df1=2,df2=N-3)
  return(1 - pf(q=f0,df1=2,df2=N-3,ncp=evalWorstNCP(mt,freqs,Amp,alpha)))
}