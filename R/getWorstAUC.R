getWorstAUC = function(tvec,Nmc,p_osc,freq,Amp,Nacro){
  acro = seq(0,2*pi,length.out=Nacro+1)
  acro = acro[1:Nacro]
  acro %>% sapply(FUN=function(phi){getAUC(tvec,Nmc,p_osc,freq,Amp,phi)}) %>% min()
}