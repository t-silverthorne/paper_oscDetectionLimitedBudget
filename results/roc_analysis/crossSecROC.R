require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(lomb)
require(pROC)
require(ggplot2)
sols = readRDS('../hiresSols.RDS')

df=readRDS('results/roc_analysis/revisedROC.RDS')

Nacro = 2^3
freqs = c(5,5.65)
acros = seq(0,2*pi,length.out=Nacro+1)
acros = acros[1:Nacro]
df=df[df$Amp==1 & df$Nmeas==40 &df$freq%in% freqs,] 

pars  = expand.grid(freq=freqs,
                    acro=acros,
                    Nmeas=c(40),
                    Amp = c(1),
                    p_osc = c(0.5),
                    fdr_method=c('none'),
                    type=c('equispaced','threshold','balanced','regu_no_cstr'))
dim(pars)

df=c(1:dim(pars)[1]) %>% lapply(function(ind){#parallel inside
  freq  = pars[ind,]$freq
  Amp   = pars[ind,]$Amp
  p_osc = pars[ind,]$p_osc
  Nmeas = pars[ind,]$Nmeas
  roc_method = pars[ind,]$roc_method
  
  mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
  mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
  mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
  mt_rnc     = sols[sols@wreg==1 & sols@drts==Inf & sols@Nmeas==Nmeas,]
  mt_opt     = mt_opt[mt_opt<Inf]
  mt_rob     = mt_rob[mt_rob<Inf]
  mt_rnc     = mt_rnc[mt_rnc<Inf]
    
  fdr_method=toString(pars[ind,]$fdr_method)
  if (pars[ind,]$type=='equispaced'){
    tvec = mt_unif
  }else if(pars[ind,]$type=='threshold'){
    tvec = mt_opt
  }else if(pars[ind,]$type=='balanced'){
    tvec = mt_rob
  }else if(pars[ind,]$type=='regu_no_cstr'){
    tvec = mt_rnc
  }else{
    stop('unrecognized type')
  }

}