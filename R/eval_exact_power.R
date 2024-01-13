#' Evaluate exact power of harmonic regression hypothesis test.
#' 
#' @description
#' Returns exact power of harmonic regression hypothesis test.
#' 
#' @param t measurement schedule
#' @param param list of oscillation parameters: 
#' \itemize{
#' \item \code{Amp}  amplitude, 
#' \item \code{freq} frequency,
#' \item \code{acro} acrophase.
#' }
#' @param alpha=.05 type I error
#' 
#' @return statistical power for given parameters
#' 
#' @author Turner Silverthorne
#' @export
eval_exact_power <- function(t,param,alpha=.05){
# return power of one-frequency cosinor model
  Amp    = param[['Amp']];
  freq   = param[['freq']];
  acro   = param[['acro']];
  N      = length(t)
  
  cvec   = Amp*cos(2*pi*freq*t-acro);
  lambda = as.numeric(t(cvec)%*%cvec)
  
  f0     = qf(p=1-alpha,df1=2,df2=N-3)
  return(1 - pf(q=f0,df1=2,df2=N-3,ncp=lambda))
}