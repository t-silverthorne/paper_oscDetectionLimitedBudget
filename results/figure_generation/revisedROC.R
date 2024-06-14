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

mc_cores = 8
Nmc      = 2e3
Nmeas    = 48

freq_vals  = seq(1,24,.5)
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
mt_rand    = runif(Nmeas)
pars       = expand.grid(freq=freq_vals,
                         Amp = c(1,1.5),
                         p_osc = c(0.1,0.5),
                         fdr_method=c('none','fdr'),
                         type=c('equispaced','threshold','balanced','random'))
pars=rbind(pars,pars,pars) # run 3 copies so you can compute sdev

df=c(1:dim(pars)[1]) %>% lapply(function(ind){
  freq  = pars[ind,]$freq
  Amp   = pars[ind,]$Amp
  p_osc = pars[ind,]$p_osc
  fdr_method=toString(pars[ind,]$fdr_method)
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
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=8,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
    return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
  }) %>% rbindlist() %>% data.frame()
    
  # FDR correction
  qval   = p.adjust(pvdf$pval,method=fdr_method)
  ostate = pvdf$state
  
  # record AUC
  roc=pROC::roc(as.numeric(ostate=='osc'),qval,direction='>')
  # record FPR and TPR at alpha=.05
  num_P   = sum(ostate=='osc') 
  num_N   = sum(ostate=='non_osc') 
  num_TP  = sum(ostate=='osc'     & qval < .05) 
  num_FP  = sum(ostate=='non_osc' & qval < .05)
  TPR     = num_TP/num_P
  FPR     = num_FP/num_N
  
  return(cbind(pars[ind,],data.frame(AUC=roc$auc,TPR=TPR,FPR=FPR)))
}) %>% rbindlist() %>% data.frame()

saveRDS(df,'results/data/roc2_hires.RDS')
# output
#df %>% ggplot(aes(x=freq,y=AUC,color=type,group=type))+
#  geom_line()+geom_point()+facet_wrap(~Amp+p_osc,ncol=4)+ylim(c(0.25,1))
#
#df %>% ggplot(aes(x=freq,y=TPR,color=type))+geom_point()+facet_wrap(~Amp+p_osc)
#df %>% ggplot(aes(x=freq,y=FPR,color=type))+geom_point()+facet_wrap(~Amp+p_osc)







