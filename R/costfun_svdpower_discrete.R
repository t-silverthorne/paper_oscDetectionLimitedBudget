costfun_svdpower_discrete = function(xinds,tau,freqs,Amp=1,alpha=.05,regL1=0){
  return(costfun_svdpower(tau[xinds],freqs,Amp,alpha,regL1))
}