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
require(matrixTests)
load_all()
sols = readRDS('results/data/MCperiodogram/hiresSols.RDS')


#sols = readRDS('../hiresSols.RDS')
#mc_cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
Nmc      = 1e4

rowCosinor <- function(theData, zts, per=24) {
  Y <- as.matrix(theData)
  
  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)
  
  betas <- solve(t(X) %*% X) %*% t(X) %*% t(Y)
  
  phases     <- atan2(betas[2,], betas[3,]) %% (2*pi)
  amplitudes <- sqrt(betas[2,]*betas[2,] + betas[3,]*betas[3,])
  
  fits <- t(X %*% betas)
  
  SStot <- rowSums((Y - rowMeans(Y))^2)
  SSres <- rowSums((fits-Y)^2)
  Rsqs  <- 1 - (SSres/SStot)
  
  SSmod <- SStot - SSres
  DFres <- ncol(theData) - 3
  DFmod <- 2
  MSres <- SSres / DFres
  MSmod <- SSmod / DFmod
  Fstatistic <- MSmod / MSres
  
  pval <- pf(Fstatistic, DFmod, DFres, lower.tail=FALSE)
  
  data.frame(phase=phases, amplitude=amplitudes, mesor=betas[1,],
             rsq=Rsqs, statistic=Fstatistic, pvalue=pval
  )
}

#high_freq=TRUE
#if (high_freq){
#  freq_vals  = seq(24,36,.05)
#  freq_vals  = freq_vals[2:length(freq_vals)]
#  fname = 'revisedROC_highfreq.RDS'
#}else{
#  freq_vals  = seq(1,24,.05)
#  fname = 'revisedROC.RDS'
#}
#pars       = expand.grid(freq=freq_vals,
#                         Nmeas=c(32,40,48),
#                         Amp = c(1,2),
#                         p_osc = c(0.5),
#                         fdr_method=c('none'),
#                         type=c('equispaced','threshold','balanced','regu_no_cstr'))

eps=1e-3
freq_vals = seq(1,36,.05)
freq_vals = freq_vals + eps*runif(length(freq_vals))
fname = 'cosinorROC.RDS'

pars       = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1),
                         p_osc = c(0.5),
                         fdr_method=c('none'),
                         stat_method=c('cosinor'),
                         type=c('equispaced','threshold','balanced','regu_no_cstr'))
 

#pars=rbind(pars,pars,pars) # run 3 copies so you can compute sdev
dim(pars)
df=c(1:dim(pars)[1]) %>% lapply(function(ind){#parallel inside
  print(ind)
  freq  = pars[ind,]$freq
  Amp   = pars[ind,]$Amp
  p_osc = pars[ind,]$p_osc
  Nmeas = pars[ind,]$Nmeas
  stat_method= pars[ind,]$stat_method
  
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
  
  # simulate data
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
  
  # simualte p-values
  if (stat_method=='cosinor'){
   pvdf = data.frame(pval=rowCosinor(Ydat,tvec,per=1/freq) %>% {.$pvalue},
                     state=state)
  }else if(stat_method=='lomb'){
    pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
      x              = Ydat[ii,]
      lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
      return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
    }) %>% rbindlist() %>% data.frame()
  }else{
    stop('unknown stat method')
  }
    
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

saveRDS(df,paste0('results/roc_analysis/',fname))
#df=readRDS('results/data/roc.RDS')
#df.sum=df %>% filter(type!='random' ) %>% group_by(freq,Nmeas,Amp,p_osc,fdr_method,type) %>% 
#  summarise(sd_AUC = sd(AUC),
#            sd_TPR = sd(TPR),
#            sd_FPR = sd(FPR),
#            AUC = mean(AUC),
#            TPR=mean(TPR),
#            FPR=mean(FPR))
#df.sum %>%  ggplot(aes(x=freq,y=AUC,color=type,group=type))+
#  geom_line()+#geom_errorbar(aes(ymin=AUC-sd_AUC,ymax=AUC+sd_AUC),data=df.sum)+
#  facet_grid(Nmeas~Amp)
#
#df.sum %>%  ggplot(aes(x=freq,y=TPR,color=type,group=type))+
#  geom_line()+#geom_errorbar(aes(ymin=AUC-sd_AUC,ymax=AUC+sd_AUC),data=df.sum)+
#  facet_grid(Nmeas~Amp)





