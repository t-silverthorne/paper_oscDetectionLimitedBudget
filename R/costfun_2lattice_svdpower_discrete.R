costfun_2lattice_svdpower_discrete = function(N1,dx1,N2,dx2,xshift2,
                                              tau,freqs,Amp=1,alpha=.05){
xinds=unique(xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2))
costfun_svdpower_discrete(xinds,tau,freqs,Amp,alpha)
}
  