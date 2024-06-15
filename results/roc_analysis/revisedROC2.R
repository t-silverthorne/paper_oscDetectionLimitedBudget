Nmc        = 5e2
Nmeas      = 48
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
mt_rand    = runif(Nmeas)

mt_opt=mt_opt[mt_opt<Inf]
mt_rob=mt_rob[mt_rob<Inf]

freqs=c(1,2,4,8,16,23,24)
p_osc = 0.5
Amp   = 1
pars=expand.grid(freq=freqs,Amp=Amp,p_osc=p_osc,type=c('equispaced','threshold','balanced','random'))
dim(pars)

rdf=c(1:dim(pars)[1]) %>% lapply(function(ind){

  freq  = pars[ind,]$freq
  Amp   = pars[ind,]$Amp
  p_osc = pars[ind,]$p_osc
  if (pars[ind,]$type=='equispaced'){
    tvec = mt_unif
  }else if(pars[ind,]$type=='threshold'){
    tvec = mt_opt
  }else if(pars[ind,]$type=='balanced'){
    tvec = mt_rob
  }else if(pars[ind,]$type=='random'){
    tvec = mt_rand
  }else{
    stop('unrecognized type')
  }  
  
  # simulate data
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
  
  # simualte p-values
  pvdf = c(1:dim(Ydat)[1]) %>% lapply(function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
    return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
  }) %>% rbindlist() %>% data.frame()
  
  # FDR correction
  qval   = p.adjust(pvdf$pval)
  ostate = pvdf$state
  
  # record AUC
  roc=pROC::roc(as.numeric(ostate=='osc'),qval,direction='>')
  
  return(cbind(pars[ind,],data.frame(sensi=roc$sensitivities,speci=roc$specificities)))
}) %>% rbindlist() %>% data.frame()


rdf %>% ggplot(aes(x=1-sensi,y=speci,group=type,color=type))+geom_line()+facet_wrap(~freq)