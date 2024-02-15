costfun_auglattice=function(N1,N2,shift2,scale2,freqs,...){
  mt=helper_auglattice_to_state(N1=N1,N2=N2,
                                shift2=shift2,scale2=scale2)
  return(costfun_svdpower(mt,freqs,...))
}