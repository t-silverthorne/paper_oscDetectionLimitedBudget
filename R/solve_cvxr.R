#'Solve power optimization using CVXR
#' 
#' @description
#' Performs DCP if \code{control$prob_formulation = 'max_elemwise'} or \code{'max_entries'}.
#' Performs hypervolume scalarization if \code{control$prob_formulation == 'hvolume'}.
#' 
#' @param control$Nmeas number of measurements to be collected
#' @param control$Nfine number of points in fine grid from which measurement times are to be chosen
#' @param control$Nfreq number of points to use in frequency discretisation of objective
#' @param control$fmin minimum candidate frequency
#' @param control$fmax maximum candidate frequency
#' @param control$PreSolve recommended \code{PreSolve=2}, see [CVXR::solve()] for details
#' @param control$MIPFocus recommended \code{PreSolve=3}, see [CVXR::solve()] for details
#' @param control$cvxr_verbose boolean for how much output is desired from CVXR
#' @param control$num_iter max number of CVXR iterations
#' @param control$MIPGapAbs stopping criterion for CVXR solver
#' @param control$time_limit max time allowed for computation
#' @param control$prob_formualtion either max_elemwise, max_entries, or hvolume, see description
#' @param Threads how many threads should CVXR use 
#'  
#' @return the result of calling [CVXR::solve] with specified user controls
#' 
#' @note
#' the \code{...} argument is only relevant for computing the cost function value
#' at the end of this method, it is not relevant to the optimization stage.
#' 
#' @note 
#' \code{time_limit} does not count the time taken to complete the pre-solve.
#' 
#' @note 
#' if iterations exceed \code{num_iter}, it is difficultto access output of the computation.
#' For this reason, we suggest setting \code{num_iter} to an unreasonably high 
#' value and using \code{time_limit} instead to determine the cut-off time for computation.
#' 
#' @seealso [opt_osc_power()] which contains a similar wrapper for CVXR
#' @author Turner Silverthorne
#' @export
solve_cvxr=function(control,Threads,...){
  start_time=Sys.time()
  weight=NULL
  Aquad             = make_quadmats(control)
  x                 = make_variable(control)
  if (control$prob_formulation=='max_elemwise'){
    prob = make_problem(x,Aquad,control)
  }else if (control$prob_formulation=='max_entries'){
    prob = make_problem_entries(x,Aquad,control) 
  }else if (control$prob_formulation=='hvolume'){
    weight=rnorm(length(Aquad)) # sample positive orthant of unit sphere
    weight=abs(weight/sqrt(sum(weight^2))) #TODO: verify
    prob = make_problem_hvolume(x,Aquad,control,weight)
  }else if(control$prob_formulation=='L1'){
    weight=rep(1,length(Aquad))/length(Aquad)
    prob = make_problem_hvolume(x,Aquad,control,weight)    
  }else{
    stop('unrecognised problem formulation')
  }
  xout = CVXR::solve(prob,verbose=control$cvxr_verbose,num_iter=control$maxit,solver="GUROBI",
                       TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,
                       Presolve=control$PreSolve,MIPFocus=control$MIPFocus,
                     Threads=Threads,NodefileStart=control$NodefileStart)
  
  tau = c(1:control$Nfine)/control$Nfine - 1/control$Nfine
  mtvalue    = tau[as.logical(xout[[1]]>1-1e-6)] # TODO: better way of catching this 
  
  freqs      = seq(from=control$fmin,to=control$fmax,length.out=control$Nfreq)
  fvalue     = costfun_svdpower(mt=mtvalue,freqs=freqs,...) 
  xindsvalue = xout[[1]]
  
  tstamp   = lubridate::now() %>% toString() %>% str_replace(' ','___')
  
  end_time =Sys.time()
  
  res_full = list(fvalue     = fvalue,
       mtvalue    = mtvalue,
       xindsvalue = xindsvalue,
       timestamp  = tstamp,
       runtime    = end_time-start_time,
       optim_raw  = xout,
       weight     = weight)
}
