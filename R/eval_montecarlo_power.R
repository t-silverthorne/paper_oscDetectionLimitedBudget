#' Use Monte Carlo method to estimate power 
#' 
#' @description
#' Performs harmonic regression on a simulated dataset of independent samples
#' and returns the portion of samples that have statistically significant p-values.
#' Useful for comparison with the exact expression for statistical power,
#' see [eval_exact_power()].
#' 
#' @param tvec vector of measurement times
#' @param param oscillation parameters
#' \itemize{
#' \item \code{Amp} amplitude of oscillation,
#' \item \code{freq} frequency of oscillation,
#' \item \code{acro} acrophase of oscillation.
#' }
#' @param Nmc number of Monte Carlo samples
#' @param alpha type I error, default value \code{alpha=.05} 
#' 
#' @return Monte Carlo estimate of power
#' @author Turner Silverthorne
#' @export
eval_montecarlo_power<-function(tvec,param,Nmc,alpha=.05){
  Amp  = param$Amp
  freq = param$freq
  acro = param$acro
  Ydat = replicate(Nmc,{Amp*cos(2*pi*freq*tvec -acro) + rnorm(length(tvec))}) %>% t
  return(rowCosinor(Ydat,tvec,per=1/freq) %>% {.$pvalue <.05} %>% mean())
}