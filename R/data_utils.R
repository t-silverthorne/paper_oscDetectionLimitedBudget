#' (Depracated) Generate simulated dataset with amplitude and freq uncertainty
#' 
#' @description
#' Generated simulated dataset with uniformly random frequency and amplitude and 
#' Gaussian white noise. The amplitude is assumed to be in units such that the noise
#' is of strength 1. In other words, the amplitude is the signal to noise ratio.
#' 
#' @param mt vector of measurement times
#' @param Nmc size of simulated dataset 
#' @param Amin minimum amplitude
#' @param Amax maximum amplitude
#' @param fmin minimum frequency
#' @param fmax maximum frequency
#' 
#' @return an [annmatrix::annmatrix()] object with rows labeled by amplitude,
#'  frequency, and acrophase. Columns are labeled by measurement time \code{mt}.
#'  
#' @author Turner Silverthorne  
#' @export
make_simulated_data=function(mt,Nmc,Amin,Amax,fmin,fmax){
  stop("depracated")
  am      = as.annmatrix(matrix(replicate(Nmc*length(mt),{0}),nrow=Nmc,ncol=length(mt)))
  am@Amp  = runif(Nmc,min=Amin,max=Amax)
  am@freq = runif(Nmc,min=fmin,max=fmax)
  am@acro = runif(Nmc,min=0,max=2*pi)

  acromat = am@acro%*%matrix(replicate(length(mt),{1}),ncol=length(mt))
  Ampmat =  am@Amp%*%matrix(replicate(length(mt),{1}),ncol=length(mt))
  epsmat  = matrix(rnorm(Nmc*length(mt)),nrow=Nmc)
  am[,]   = Ampmat*cos(2*pi*am@freq%*%t(mt)-acromat)
  am[,]   = am[,]+epsmat 
  am$time = mt
  return(am)
}

