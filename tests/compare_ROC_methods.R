require(devtools)
load_all()
require(pROC)
require(ggplot2)

Nmc        = 1e4
Nmeas      = 20
tvec       = c(1:Nmeas)/Nmeas-1/Nmeas
alpha_vals = seq(0,1,.01)
sim        = simLombScarglePvals(Nmc,tvec,p_osc=0.5)
roc1       = getLombScargleROC(alpha_vals,Nmc,tvec,p_osc=0.5) %>% filter(p_method=='std' & fdr_method=='none')

pvec   = sim$pvdf %>% filter(p_method=='std') %>% {.$pval}
bstate = sim$pvdf %>% filter(p_method=='std') %>% {.$state} %>% 
  sapply(function(x){ifelse(x=='osc',1,0)}) %>% matrix()

head(bstate)
head(pvec)
roc2=pROC::roc(factor(bstate),pvec,plot=T)


df1=data.frame(x=roc1$FPR,
           y=roc1$TPR,
           method='roc1')

df2=data.frame(x=1-roc2$specificities,
           y=roc2$sensitivities,
           method='roc2')

rbind(df1,df2) %>% ggplot(aes(x=x,y=y,group=method,color=method))+geom_line()