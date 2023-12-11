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
  Tinit        = 100
  Tnow         = Tinit 
  # initialize state 
  if (opts$lattice_cstr=='none'){
    x = rep(0,opts$Nfine)
    x[sample(c(1:opts$Nfine),opts$Nmeas,T)]=1
  } else if (opts$lattice_cstr=='sa_lattice'){
    x = generate_random_lattice(opts,x=NULL)
  } else { 
    stop('Unrecognized lattice constraint. Did you forget to change from CVXR options?')
  }
  ii=1
  Sx=sa_cfunpwr(x,Aquad,opts)
  while (ii<opts$num_iter+1){
    if (ii %% floor(opts$num_iter/20) == 0){
      print(paste0('Completed: ', toString(100*ii/opts$num_iter),' perc of SA run. Temp: ', toString(Tnow)))
      #Tnow  = Tinit*.65
      #Tinit = Tnow
    }
    y = sa_propfunction(x,opts)
    
    Sy=sa_cfunpwr(y,Aquad,opts)      # eval cost fun of candidate state
    alpha = min(exp(-(Sy-Sx)/Tnow),1)         # acceptance prob
    
    Tnow = Tnow*cooling_rate
    if (rands[ii]<alpha){
      x  = y
      Sx = Sy
    }
    ii=ii+1
  }
  return(Sx) 
}
