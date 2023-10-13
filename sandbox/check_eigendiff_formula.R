# verify that eigenvalue formula matches finite difference
library(stats)
library(ggplot2)
library(dplyr)
library(expm)
library(patchwork)
library(pracma)
library(data.table)
source('utils/powerUtils.R')

N=12
param=list(Amp=1.5,
           freq=4.4,
           acro=1,
           method='min')
t=(0:N)/N
t=t[1:(length(t)-1)]
t0=t

fmin = 1 
fmax = 24
df   = .1 

freqlist=seq(from=fmin,to=fmax,by=1e-2)
df = freqlist %>% lapply(function(freq){
  # exact formula
  param$freq=freq
  dLam = diffLambdaMin(t,param)
  lam  = getMinEig(t,param)
  
  # finite difference
  dfreqfd = 1e-5
  param$freq=freq+dfreqfd
  lam2   = getMinEig(t,param)
  dLamFD = (lam2-lam)/dfreqfd
  return(rbind(data.frame(freq=freq,eig=lam,dLam=dLam,method='exact'),
               data.frame(freq=freq,eig=lam,dLam=dLamFD,method='FD')))
}) %>% rbindlist() 


df %>% ggplot(aes(x=freq,y=dLam,group=method,color=method))+geom_line()+facet_wrap(~method,nrow=2)