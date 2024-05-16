require(devtools)
require(ggplot2)
devtools::load_all()
Nmeas = 48
tvec  = c(1:Nmeas)/Nmeas-1/Nmeas
df    = getLombScargleROC(seq(0,1,.05),1e3,tvec,freq=1)
df %>% ggplot(aes(x=FPR,y=TPR))+geom_line()+facet_grid(p_method~fdr_method)
