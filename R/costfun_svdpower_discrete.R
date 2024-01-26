costfun_svdpower_discrete = function(xinds,Nfine,freqs,Amp=1,alpha=.05){
  tau=seq(from=0,to=1,length.out=Nfine+1)
  tau=tau[1:(length(tau)-1)]
  tau=tau[xinds]
  
  return(costfun_svdpower(tau,freqs,Amp,alpha))
}