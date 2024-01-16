#' Simulated annealing routine for lattice constraints in power maximization
#' 
#' @description
#' Use simulated annealing to generate candidate solutions to power maximization problem
#'
#' @param Aquad the matrix
#' @param opts options
#' 
#' @return simulated annealing result
#'
#' @author Turner Silverthorne
#' @export
run_sa_power=function(Aquad,opts){
  if (!is.null(attributes(class(Aquad)))){
    if (attr(class(Aquad),'package')=='CVXR'){
      stop('You are running simulated annealing, Aquad should not be a CVXR constant')
    }
  }
  if(any(class(Aquad)=='list')){
    if (!is.null(attributes(class(Aquad[[1]])))){
      if (attr(class(Aquad[[1]]),'package')=='CVXR'){
        stop('You are running simulated annealing, first entry of Aquad should not be a CVXR constant')
      }
    }
  }
 
  #TODO: make this user determined
  epochs       = F 
  rands        = runif(opts$num_iter)
  cooling_rate = .95
  Tinit        = 100
  Tnow         = Tinit 
  
  # initialize state 
  if (opts$lattice_cstr %in% c('none','sa_lattice')){
    x = sa_propfunction(opts,x=NULL) # random initialization 
  }else{
    stop('Unrecognized lattice constraint. Did you forget to change from CVXR options?')  
  }
  
  
  # run annealing   
  ii=1
  Sx=sa_cfunpwr(x,Aquad,opts)
  while (ii<opts$num_iter+1){
    if (ii %% ceiling(opts$num_iter/20) == 0){
      if (opts$verbose){
        print(paste0('Completed: ',
                     toString(round(100*ii/opts$num_iter,2)),
                     ' perc of SA run. Temp: ',
                     toString(round(Tnow,2)),
                     ' fval: ',toString(round(Sx,4))))
      }
      if (epochs){
        Tnow  = Tinit*.65
        Tinit = Tnow
      }
    }
    y = sa_propfunction(opts,x)
    
    Sy=sa_cfunpwr(y,Aquad,opts)      # eval cost fun of candidate state
    alpha = min(exp(-(Sy-Sx)/Tnow),1)         # acceptance prob
    
    Tnow = Tnow*cooling_rate
    if (rands[ii]<alpha){
      x  = y
      Sx = Sy
    }
    ii=ii+1
  }
  return(list(xval=x,cfunval=Sx))
}
