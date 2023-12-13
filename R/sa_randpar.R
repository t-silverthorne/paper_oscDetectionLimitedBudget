#' Generate uniformly random partition of integer n
#' 
#' @description
#' This function returns a uniformly random partition of the input integer \code{n}. The algorithm here
#' is from Combinatorial Algorithms by Nijenhuis and Wilf (1975). 
#' 
#' It is a useful for constructing the transition function for running simulated annealing 
#' with lattice constraints.
#' 
#' @param n non-negative integer which one wishes to partition
#' 
#' @return a uniformly random partition of \code{n}. 
#' 
#' @author Turner Silverthorne
sa_randpar=function(n){
  P      = list()
  m      = n
  while (m>0){
    emat   = sa_eulermat(m) 
    states = lapply(emat,function(x){x[1:2]})
    probs  = unlist(lapply(emat,function(x){x[3]}))
    state = sample(states,1,FALSE,probs)
    dd=state[[1]][1]
    jj=state[[1]][2]
    P=append(P,replicate(jj,{dd}))
    m=m-jj*dd
  }
  return(unlist(P))
}