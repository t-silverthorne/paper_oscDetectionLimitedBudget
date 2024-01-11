#' Power maximization using disciplined convex programming
#' 
#' @description
#' Use [CVXR] with [GUROBI] backend to construct optimal sampling strategies for
#' fixed sampling budget and frequency uncertainty.
#' 
#' For lattice constraints, it is usually much more efficient to use the simulated
#' annealing solver [run_sa_power()].
#' 
#' @param prob a convex programming problem constructed using [make_problem()]
#' @param opts optimization options, for instance from [make_default_opts()]
#' 
#' @return disciplined convex programming result
#' 
#' @author Turner Silverthorne
run_cvxr_power=function(prob,opts){
  #start    = Sys.time()
  result   = CVXR::solve(prob,verbose=opts$verbose,num_iter=opts$num_iter,
                   MIPGapAbs=opts$MIPGapAbs)
  #end      = Sys.time()
  return(result)
}

