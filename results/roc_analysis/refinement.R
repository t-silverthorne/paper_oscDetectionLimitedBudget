require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(ggplot2)
load_all()
mc.cores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
rowCosinor <- function(theData, zts, per=24) {
  Y <- as.matrix(theData)
  
  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)
  
  betas <- solve(t(X) %*% X,t(X) %*% t(Y))
  
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


mt = c(1:16)/16-1/16
mt = .25*mt

freqs = c(3.9,3.5) 
Nmc   = 1e4
nreps  = c(1e5,1e4,1e3,1e2,1e1)

Nmeas=length(mt)
Amp=1.2
acro=0
pars=expand.grid(freq=freqs,nrep=nreps)
df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=12,function(ind){
  freq=pars[ind,]$freq
  nrep=pars[ind,]$nrep
  
  param=list(freq=freq,acro=acro,Amp=Amp)
  pwr_exact=evalExactPower(mt,param) 
 
  pwr_est=0 
  for (ii in c(1:nrep)){
    Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)+
      t(replicate(Nmc,Amp*cos(2*pi*mt*freq-acro)))
    pwr_est = pwr_est+(rowCosinor(Ydat,mt,1/freq) %>% {.$pvalue<.05} %>% mean())/nrep
  }
  
 return(cbind(pars[ind,],data.frame(err=abs(pwr_est-pwr_exact)/pwr_exact) ))
}) %>% rbindlist() %>% data.frame()
saveRDS(df,'results/roc_analysis/refinement.RDS')
#df %>% ggplot(aes(x=nrep*Nmc,y=err,group=freq,color=freq))+geom_line()+
#  scale_x_continuous(trans='log10')+
#  scale_y_continuous(trans='log10')