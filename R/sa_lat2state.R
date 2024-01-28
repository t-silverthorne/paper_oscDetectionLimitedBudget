#' Convert lattice structure to binary vector representing measurements
#' 
#' @description
#' When using lattice constraints for simulated annealing, the state vector is
#' represented by a list of lattices. This function takes that list and constructs
#' a binary vector of length \code{opts$Nfine} representing the actual state of
#' the system.
#' 
#' This function is mostly used as a helper function for evaluating the simulated 
#' annealing cost function found in [sa_cfunpwr()].
#' 
#' @param x current state of system
#' @param Nfine number of points in time discretization
#' 
#' @author Turner Silverthorne
#' @export
sa_lat2state=function(x,Nfine){
  y = rep(0,Nfine)
  for (ii in 1:length(x)){
    xloc = x[[ii]]
    y = y + xloc
  }
  return(y)
}