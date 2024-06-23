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