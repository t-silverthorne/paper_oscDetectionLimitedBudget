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

alpha_vals = seq(0,1,.05)
freq_vals  = c(1,2,4,8,16,24) 
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
pars       = expand.grid(freq=freq_vals,
                         Amp = c(1,1.5),
                         p_osc = c(0.1,0.5),
                         type=c('equispaced','threshold','balanced'))

dfROC = c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ind){
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
  
  TPR_loc = getLombScargleROC(alpha_vals = alpha_vals,
                        Nmc   = Nmc,
                        tvec  = tvec_loc,
                        Amp   = Amp,
                        freq  = freq,
                        p_osc = p_oscloc,
                        mc_cores=1)
  return(cbind(pars[ind,],TPR_loc))
}) %>% rbindlist() %>% data.frame()


dfROC %>% filter(fdr_method=='BH',p_method=='press') %>% 
  ggplot(aes(x=FPR,y=TPR,group=type,color=type))+geom_line()+
  scale_x_continuous(limits = c(0,1))+
  facet_grid(Amp+fdr_method+p_osc+p_method~freq)