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
  if (is.null(opts)){  # unpack options
    lattice_cstr='none'
  }else if(opts$lattice_cstr %in% c('none','sa_lattice')){
    lattice_cstr=opts$lattice_cstr
  }else{
    stop('Unrecognized lattice constraint')
  }
  
  if (is.null(parent_size)){  # unpack parent_size
    parent_size = 0
  }
  
  if (is.null(opts) | lattice_cstr=='none'){
    if (parent_size>0){
      stop('Provided size of parent lattice without enforcing lattice constraints. This 
           makes no sense because parent lattice size is only relevant to partition 
           generation when the lattice constraints are being enforced.')
    }
  } 
  
  if (!is.null(opts)){
    if (lattice_cstr=='none'){
      warning('You provided options to sa_randpar but they will be ignored since
              your opts$lattice_cstr=none. Did you mean to specify lattice constraints?')
    }else if(lattice_cstr=='sa_lattice'){ # check at least one partition can satisfy constraints 
      if(n < opts$min_lat){
        stop('Number of points to be partitioned is less than min lattice size') 
      }
      if(n+parent_size < opts$min_active_lats){
        stop('No feasible partition: you must decrease min_active_lats') 
      }
      if(1+parent_size > opts$max_active_lats){
        stop('No feasible partition: you must increase max_active_lats') 
      }
      if (opts$max_lat*opts$max_active_lats < n){
        stop('Problem is infeasible: increase number of allowed lattices and/or 
        max number of allowed points in lattice')
      }
      if (opts$min_lat*opts$min_active_lats > n){
        stop('Problem is infeasible: decrease min number of allowed lattices and/or 
          min number of allowed points in lattice')
      }
    }
  }
  
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
    if (lattice_cstr=='none'){
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
