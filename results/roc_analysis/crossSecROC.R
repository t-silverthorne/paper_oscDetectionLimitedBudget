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
load_all()
sols = readRDS('results/data/MCperiodogram/hiresSols.RDS')

df=readRDS('results/roc_analysis/revisedROC.RDS')

mc_cores = 12
Nmc      = 5e2
Nacro    = 2^5
freqs    = c(5,5.65)
acros    = seq(0,2*pi,length.out=Nacro+1)
acros    = acros[1:Nacro]
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
  acro  = pars[ind,]$acro
  
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
  return(cbind(pars[ind,],data.frame(AUC=getAUC(tvec  = tvec,
                                                Nmc   = Nmc,
                                                p_osc = p_osc,
                                                freq  = freq,
                                                Amp   = Amp,
                                                acro  = acro))))
}) %>% rbindlist() %>% data.frame()
saveRDS('results/roc_analysis/crossSecROC.RDS')

#df=readRDS('results/roc_analysis/crossSecROC.RDS')
#df %>% ggplot(aes(x=acro,y=AUC,group=freq,color=freq))+geom_line()+facet_wrap(~type)
