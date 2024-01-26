#test_that("function evaluation", {
#  library(annmatrix)
#  library(dplyr)
#  library(data.table)
#  mt   = c(0:15)/15
#  mt   = mt[1:(length(mt)-1)]
#  Nmc  = 5e2
#  fmin = 1
#  fmax = 24 
#  Amin = 0.5 
#  Amax = 5 
#  Namp = 10 
#  true_freq_vals  = c(1,2,12) %>% as.list()
#  Ampvals         = seq(from=Amin,to=Amax,length.out=Namp) %>% as.list()
#  Nfreq_regr_vals = c(10,15,20) %>% as.list()
#  start=Sys.time()
#  uu=benchmark_fdr(mt,Nmc,fmin,fmax,Ampvals,Nfreq_regr_vals,'fdr')
#  end=Sys.time()
#  end-start
#  library(ggplot2)
#  library(patchwork)
#  uu %>% head()
#  
#  #p1=uu %>% ggplot(aes(x=Amp,y=pdetect_p,color=as.factor(true_freq),group=true_freq))+
#  #  geom_line()+facet_wrap(~Nfreq_regr,nrow=3) +theme(legend.position='none')
#  #p2=uu %>% ggplot(aes(x=Amp,y=pdetect_q,color=as.factor(true_freq),group=true_freq))+
#  #  geom_line()+facet_wrap(~Nfreq_regr,nrow=3) +theme(legend.position='none')
#  #
#  #p1+p2
#})
