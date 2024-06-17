sim_data = function(tvec){
  mc_cores=10
  fdr_method='none'
  Nmc   = 1e4
  Nmeas = length(tvec)
  p_osc = 0.5
  Amp   = 2
  freq  = 1 
  
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
  return(list(qval=qval,ostate=ostate))  
}

Nmeas=48
tvec_unif  = c(1:Nmeas)/Nmeas - 1/Nmeas
tvec       = sols[sols@wreg==1 & sols@drts==Inf & sols@Nmeas==Nmeas,]
tvec=tvec[tvec<Inf]

out1=sim_data(tvec_unif)
out2=sim_data(tvec)
# record AUC
roc1=pROC::roc(as.numeric(out1$ostate=='osc'),out1$qval,direction='>',plot=T)
roc2=pROC::roc(as.numeric(out2$ostate=='osc'),out2$qval,direction='>',plot=T)
roc1
roc2
?roc.test()