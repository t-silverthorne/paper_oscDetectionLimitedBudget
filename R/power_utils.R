#' Evaluate exact power of harmonic regression hypothesis test.
#' 
#' @description
#' Returns exact power of harmonic regression hypothesis test. 
#' 
#' 
#' @param t measurement schedule
#' @param param$Amp amplitude of signal 
#' @param param$freq frequency fo signal
#' @param param$acro phase of signal in radians
#' @param alpha type I error, by default \code{alpha=.05}
#' 
#' @note
#' Assumes the noise has mean zero and unit standard deviation.
#' In these units, amplitude can be interpreted as the signal to noise ratio.
#' 
#' @return statistical power 
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

#'Use Monte Carlo method to estimate power 
#' 
#' @description
#' Performs harmonic regression on a simulated dataset of independent samples
#' and returns the portion of samples that have statistically significant p-values.
#' Useful for comparison with the exact expression for statistical power,
#' see [eval_exact_power].
#' 
#' @param tvec vector of measurement times
#' @param param$Amp amplitude of signal 
#' @param param$freq frequency fo signal
#' @param param$acro phase of signal in radians
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

#' Calculate FDR-corrected power for a range of phases and frequencies
#' 
#' @description
#' Given a measurement schedule, true frequency, and amplitude, the FDR-corrected power
#' is estimated for each acrophase and the minimum power is returned.
#' 
#' @note
#' Multiple test correction is performed only with respect to the multiple frequencies and
#' not with respect to the multiple phases.
#' 
#' @param mt vector of measurement times
#' @param Amp amplitude of signal
#' @param f0 true frequency of signal
#' @param Nacro number of acrophases to use in simulating the power
#' @param Nfreq number of candidate frequencies to try in harmonic regression
#' @param fmin lowest candidate frequency
#' @param fmax highest candidate frequency
#' @param fdr_method passed to the [p.adjust] function, this determines what type
#' of multiple test correction should be performed
#' 
#' @return the worst-case FDR-corrected power
#' @author Turner Silverthorne
#' @export
benchmark_worst_fdr=function(mt,Amp,f0,Nmc,Nacro,Nfreq,fmin,fmax,
                             fdr_method='BH'){
    freqs = seq(from=fmin,to=fmax,length.out=Nfreq)
    acros = seq(from=0,to=2*pi,length.out=Nacro)
    
    acros %>%sapply(function(acro){
      X=t(replicate(Nmc,{Amp*cos(2*pi*f0*mt-acro)}))+
        matrix(rnorm(Nmc*length(mt)),nrow=Nmc)
     
      Lfit = freqs %>%as.list() %>% lapply(function(freq){
        rowCosinor(X,zts=mt,per=1/freq) %>% {data.frame(pvalue=t(.$pvalue))}})
     
      pmat=Lfit%>% rbindlist() %>% as.matrix() %>% t()
      qmat=matrix(replicate(Nmc*Nfreq,{NaN}),nrow=Nmc)
      for (ii in c(1:Nmc)){
        qmat[ii,]=p.adjust(pmat[ii,],fdr_method)
      }
      is_osc = apply(qmat,1,function(x){any(x<.05)}) 
      if (length(is_osc)==Nmc){
        return(mean(is_osc))
      }else{
        stop('inconsistent dimension, did you transpose something?')
      }
    }) %>% min()
}

