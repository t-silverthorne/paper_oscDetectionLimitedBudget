require(dplyr)
require(ggplot2)
require(data.table)
require(parallel)
Nrep = 1e3
take_mean=F
u=c(1:Nrep) %>% lapply(function(ii){Nmc=10
  pvec  = runif(Nmc) # TODO: harmonic regression
  BHvec = p.adjust(pvec,method='BH')
  BYvec = p.adjust(pvec,method='BY')
  if (take_mean){
    BHvec = mean(BHvec)
    BYvec = mean(BYvec)
  }
  return(data.frame(BHval=BHvec,BYval=BYvec))
}) %>% rbindlist()

theme_set(theme_classic())
u %>% 
  ggplot(aes(x=BHval,y=BYval)) +geom_point(size=.5,alpha=.5)+
  geom_abline(slope=1,intercept=0,color='grey')+xlim(c(0,1))+ylim(c(0,1))