require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(dplyr)
require(lomb)
require(pROC)
require(ggplot2)

mc_cores   = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
sols = readRDS('../hiresSols.RDS')
getAUC = function(tvec,Nmc,p_osc,freq,Amp,acro){
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',] +
    t(replicate(sum(state=='osc'),Amp*cos(2*pi*tvec*freq-acro)))
  
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
    return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
  }) %>% rbindlist() %>% data.frame()
  
  # FDR correction
  qval   = pvdf$pval
  ostate = pvdf$state
  
  roc=pROC::roc(as.numeric(ostate=='osc'),qval,direction='>',quiet=T)
  return(roc$auc) 
}
# get worstAUC
getWorstAUC = function(tvec,Nmc,p_osc,freq,Amp,Nacro){
  acro = seq(0,2*pi,length.out=Nacro+1)
  acro = acro[1:Nacro]
  acro %>% sapply(FUN=function(phi){getAUC(tvec,Nmc,p_osc,freq,Amp,phi)}) %>% min()
}
freq_vals  = seq(1,24,.25) # todo make this a seq
Nmc        = 1e3
Nacro      = 2^5
# compare designs
pars       = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1,2),
                         p_osc = c(0.5),
                         type=c('equispaced','threshold','balanced','regu_no_cstr'))
dim(pars)

#df=c(1:dim(pars)[1]) %>% lapply(function(ind){#parallel inside
df=c(1:128) %>% lapply(function(ind){#parallel inside
  freq  = pars[ind,]$freq
  Amp   = pars[ind,]$Amp
  p_osc = pars[ind,]$p_osc
  Nmeas = pars[ind,]$Nmeas

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
  return(cbind(pars[ind,],data.frame(wAUC=getWorstAUC(tvec=tvec,
                                                      Nmc=Nmc,
                                                      p_osc=p_osc,
                                                      freq=freq,
                                                      Amp=Amp,
                                                      Nacro=Nacro))))
}) %>% rbindlist() %>% data.frame()
saveRDS(df,'results/roc_analysis/worst_case_auc.RDS')
