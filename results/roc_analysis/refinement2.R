require(matrixTests)
require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(ggplot2)
load_all()
Nmeas = 10
mt=c(1:Nmeas)/Nmeas-1/Nmeas
#mt=mt*.25
acro=0
Amp=1
freq=0.5

nrep    = 1
Nmc     = 1e3
acro = seq(0,2*pi,length.out=2^4+1)

acro = acro[c(1:2^4)]


pwr_est=c()
pwr_exact=c()
for (ii in c(1:2^4)){
  param=list(Amp=Amp,freq=freq,acro=acro[ii])
  pwr_est[ii]   = evalMonteCarloPower(mt,param,Nmc)
  pwr_exact[ii] = evalExactPower(mt,param)
}
df1=data.frame(acro=acro,pwr=pwr_est,type='est')
df2=data.frame(acro=acro,pwr=pwr_exact,type='exact')
rbind(df1,df2) %>% ggplot(aes(x=acro,y=pwr,group=type,color=type))+geom_line()+
  ylim(c(0,1))+geom_vline(xintercept = pi/2)


