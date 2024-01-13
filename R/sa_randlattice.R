#' Random lattice for simulated annealing
#'
#' @description
#' Generates a random binary lattice for the simulated annealing transition function.
#' See [run_sa_power()] for details. When there are lattice constraints, a more
#' complicated lattice proposal function is called in [sa_propfunction()].
#'  
#' @param n number of measurements
#' @param opts$Nfine coarseness of binary lattice see [make_default_opts()]
#' 
#' @author Turner Silverthorne
#' @export
sa_randlattice=function(n,opts){
  Nfine = opts$Nfine
  x     = replicate(Nfine,0)
 
  dx    = sample(1:floor(Nfine/n),1) #TODO: clarify assumption that entire lattice must fit in simulation 
  shift = sample(1:dx,1) #TODO: check that you dont need -1
  for (ii in c(1:n)){
    x[1+(shift+dx*ii) %% Nfine]=1
  }
  return(x)
}
