costfun_svdpower_discrete = function(xinds,tau,freqs,Amp=1,alpha=.05){
  return(costfun_svdpower(tau[xinds],freqs,Amp,alpha))
}