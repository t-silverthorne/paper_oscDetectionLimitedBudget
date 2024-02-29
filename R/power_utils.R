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


#' Perform multiple test correction on rows or columns of a p-value matrix
#' 
#' @description
#' User chooses if multiple test correction should be performed row-wise or column-wise.
#' Depending on this choice, rows or columns are treated as the p-values from
#' independent hypothesis tests and are corrected using the [stats::p.adjust()] function.
#' 
#' @param pdat matrix of pvalues
#' @param dim dimension along which adjustment should be performed
#' \itemize{
#' \item if \code{dim==1}, each row is treated as a collection of hypothesis tests for which
#'  a set of q-values should be computed
#'  \item if \code{dim==2}, each column is treated as a collection of hypothesis tests for which
#'  a set of q-values should be computed
#' }
#' @param pmethod adjustment method, passed to [stats::p.adjust()]
#' 
#' @return matrix of q-values
#' @author Turner Silverthorne
#' @export
matrix_1d_padjust=function(pdat,dim,pmethod){
  stop('depracated')
  qdat=NaN*pdat
  
  #TODO: could write more concisely using apply() but this gave cryptic rann errors
  if (dim==1){ 
    for (ii in c(1:dim(pdat)[1])){
      qdat[ii,] = p.adjust(pdat[ii,],method=pmethod) 
    }
  }
  
  if (dim==2){
    for (ii in c(1:dim(pdat)[2])){
      qdat[,ii] = p.adjust(pdat[,ii],method=pmethod) 
    }
  }
  return(qdat)
}