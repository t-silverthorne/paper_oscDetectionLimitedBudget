costfun_nlattice_svdpower_discrete=function(xlat_list,Nfine,freqs,Amp=1,alpha=.05){
  xinds = sa_lat2state(xlat_list,Nfine)
  costfun_svdpower_discrete(xinds,Nfine,freqs,Amp,alpha)
}