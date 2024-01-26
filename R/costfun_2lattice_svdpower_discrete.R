costfun_2lattice_svdpower_discrete = function(N1,dx1,N2,dx2,xshift2,
                                              Nfine,freqs,Amp=1,alpha=.05){

lat1  = dx1*c(0:(N1-1)) 
lat2  = xshift2 + dx2*c(0:(N2-1))
xinds = c(lat1,lat2)

costfun_svdpower_discrete(xinds,Nfine,freqs,Amp,alpha)
}
  