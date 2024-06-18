require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(dplyr)
require(annmatrix)
require(lomb)
require(pROC)
require(ggplot2)
load_all()
mc_cores=12
fdr_method='fdr'
Nmc   = 1e4
Nmeas = 48 
p_osc = 0.1
Amp   = 1 
freq  = 24
tvec  = c(1:Nmeas)/Nmeas - 1/Nmeas
#tvec   = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
tvec=tvec[tvec<Inf]
Nmeas               = length(tvec)
Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
N_osc               = sum(state=='osc')
Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))

# simualte p-values
pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
  x              = Ydat[ii,]
  lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
  return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
}) %>% rbindlist() %>% data.frame()
  
# FDR correction
qval   = p.adjust(pvdf$pval,method=fdr_method)
ostate = pvdf$state

# record AUC
pROC::roc(as.numeric(ostate=='osc'),qval,direction='>',thresholds=c(0),plot=T)
