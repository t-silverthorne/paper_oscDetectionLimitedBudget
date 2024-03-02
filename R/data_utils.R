#' (Depracated) Generate simulated dataset with amplitude and freq uncertainty
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

