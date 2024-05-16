require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(lomb)
load_all()
sols = readRDS('results/data/MCperiodogram/hiresSols.RDS')
mc_cores = 16
# redo TPR estimates
Nmc      = 1e3
Nmeas    = 48
p_osc    = .1

alpha_vals = seq(0,.25,.01)
freq_vals  = c(1,2,4,8,16,24) 
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
pars       = expand.grid(alpha_val=alpha_vals,
                         freq=freq_vals,
                         Amp = c(1,1.5,2),
                         p_osc = c(0.1,0.2,0.5),
                         type=c('equispaced','threshold','balanced'))

dfTPR = c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ind){
  alpha_val = pars[ind,]$alpha_val
  freq      = pars[ind,]$freq
  Amp       = pars[ind,]$Amp
  p_oscloc  = pars[ind,]$p_osc
  if (pars[ind,]$type=='equispaced'){
    tvec_loc = mt_unif
  }else if(pars[ind,]$type=='threshold'){
    tvec_loc = mt_opt
  }else if(pars[ind,]$type=='balanced'){
    tvec_loc = mt_rob
  }else{
    stop('unrecognized type')
  }
  
  TPR_loc = estimateTPR(alpha = alpha_val,
              Nmc   = Nmc,
              tvec  = tvec_loc,
              freq  = freq,
              Amp   = Amp,
              p_osc = p_oscloc,
              mc_cores=1)
  par_ind = rbindlist(replicate(dim(TPR_loc)[1],pars[ind,],F))
  return(cbind(par_ind,TPR_loc))
}
) %>%  rbindlist() %>% data.frame()
saveRDS(dfTPR,'results/data/MCperiodogram/dfTPR_large.RDS')