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
#' @param opts options to specify restrictions on partition, will be enforced using rejection sampling
#' @return a uniformly random partition of \code{n}. 
#' 
#' @author Turner Silverthorne
sa_randpar=function(n,opts=NULL,parent_size=NULL){
  accept_partit = F
  
  while (!accept_partit){
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
    if (is.null(opts)){
      accept_partit = T 
    }else{
      c1 = min(unlist(P))>=opts$min_lat 
      c2 = max(unlist(P))<=opts$max_lat 
      
      if (is.null(parent_size)){
        c3=length(P) >= opts$min_active_lats
        c4=length(P) <= opts$max_active_lats
      }else{
        c3=length(P)+parent_size >= opts$min_active_lats
        c4=length(P)+parent_size <= opts$max_active_lats
      }
      accept_partit = c1*c2*c3*c4
    }
  }
  return(unlist(P))
}