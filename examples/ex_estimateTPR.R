require(ggplot2)
require(dplyr)
require(devtools)
devtools::load_all()
mc_cores   = 12   # number of threads for parallel
Nmc        = 5000 # number of MC samples (<500 in parallel runs in about 1 minute)

# parameters for simulated data
p_osc      = .1   # proportion of oscillating genes
Nmeas      = 48   # number of measurements
alpha_vals = seq(0,1,.1)
freq_vals  = c(1,24) 
pars       = expand.grid(alpha=alpha_vals,
                         freq=freq_vals,
                         Amp = c(1.5),
                         type=c('uniform','optimal','robust'))

# measurement schedules
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
source('results/data/extract_data.R')
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]


df = c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
  x = pars[ii,]
  type_loc = x[['type']]
  if (type_loc =='uniform'){
    tvec_loc = mt_unif
  }else if(type_loc=='optimal'){
    tvec_loc = mt_opt
  }else if(type_loc=='robust'){
    tvec_loc = mt_rob
  }else{
    stop('unrecognized type')
  }
  
  TPRdf = estimateTPR(alpha = x[['alpha']],
              Nmc   = Nmc,
              tvec  = tvec_loc,
              freq  = x[['freq']],
              Amp   = x[['Amp']],
              p_osc = p_osc,
              mc_cores=1)
  
  DFinput = pars[ii,] 
  
  return(cbind(DFinput,TPRdf))
}) %>% rbindlist() %>% data.frame()


df %>% filter(p_method =='std') %>% 
  ggplot(aes(x=alpha,y=TPR,group=type,color=type))+
  geom_line()+
  facet_grid(fdr_method~freq)
