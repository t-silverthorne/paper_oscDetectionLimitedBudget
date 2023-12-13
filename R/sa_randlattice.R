sa_randlattice=function(n,opts){
  x     = replicate(Nfine,0)
 
  dx    = sample(1:floor(Nfine/n),1) #TODO: clarify assumption that entire lattice must fit in simulation 
  shift = sample(1:dx,1) #TODO: check that you dont need -1
  for (ii in c(1:n)){
    x[1+(shift+dx*ii) %% Nfine]=1
  }
  return(x)
}
