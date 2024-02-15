tfun_auglattice=function(N1,N2,shift2,scale2,freqs,control){
  if (N1<control$N1min){
    stop('too few points in lattice1')
  }
  if (N1>control$N1max){
    stop('too many points in lattice1')
  }
  # update number of points in each lattice
  N1cands=c(N1,min(N1+1,control$N1max),max(N1-1,control$N1min))
  N1new = sample(N1cands,1) 
  N2new = N2+N1-N1new
  
  shift2new = runif(1,0,1/N1new)
  scale2new = runif(1,0,N2new/(N2new-1)-shift2new)
  return(c(N1new,N2new,shift2new,scale2new))
}