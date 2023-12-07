#' Simulated annealing routine for lattice constraints in power maximization
#' 
#' @description
#' Use simulated annealing to generate candidate solutions to power maximization problem
#'
#' @param Aquad the matrix
#' @param opts options
#' @param Lmat constraints
#' 
#' @return simulated annealing result
#'
#' @author Turner Silverthorne

run_sa_power=function(Aquad,opts,Lmat=NULL){
  if (!is.null(Lmat)){
    if (attr(class(Lmat),'package')=='CVXR'){
      warning('You are running simulated annealing, Lmat should not be a CVXR constant')
    }
  }
  
  if (!is.null(attributes(class(Aquad))) | !is.null(attributes(class(Aquad)) )){
    if (attr(class(Aquad),'package')=='CVXR'|attr(class(Aquad[[1]]),'package')=='CVXR' ){
      warning('You are running simulated annealing, Aquad should not be a CVXR constant')
    }
  }
  
  rands        = runif(opts$num_iter)
  cooling_rate = .9999
  #Tinit        = 100
  Tnow         = Tinit 
  x = rep(0,opts$Nfine)
  # initialize state 
  if (opts$lattice_cstr=='none'){
    x[sample(c(1:opts$Nfine),opts$Nmeas,T)]=1
  } else if (opts$lattice_cstr=='lineq' | opts$lattice_cstr=='cfun'){
    stop('Simulated annealing with lattice constraints has not been implemented yet')
  }
  ii=1
  Sx=cfun_sim_anneal_pwr(x,Aquad,opts)
  while (ii<opts$num_iter+1){
    if (ii %% floor(opts$num_iter/20) == 0){
      print(paste0('Completed: ', toString(100*ii/opts$num_iter),' perc of SA run. Temp: ', toString(Tnow)))
      #Tnow  = Tinit*.65
      #Tinit = Tnow
    }
    y = x
    
    swap_on     = sample(which(x==0),1)       # proposed state differs from x at one index
    swap_off    = sample(which(x>0),1)
    y[swap_on]  = 1
    y[swap_off] = 0
    
    Sy=cfun_sim_anneal_pwr(y,Aquad,opts)      # eval cost fun of candidate state
    alpha = min(exp(-(Sy-Sx)/Tnow),1)         # acceptance prob
    
    Tnow = Tnow*cooling_rate
    if (rands[ii]<alpha){
      x=y
    }
    Sx=cfun_sim_anneal_pwr(x,Aquad,opts)
    ii=ii+1
  }
  return(Sx) 
}
