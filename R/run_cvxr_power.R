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
#' \itemize{
#' \item \code{opts$verbose} boolean for if you want solver to display progress 
#' \item \code{opts$time_limit} time limit for gurobi solver in units of seconds
#' \item \code{opts$MIPGapAbs} threshold for terminanting solver, measures gap between
#' convex relaxation of problem and the actual integer constrained programming problem
#' }
#'
#' @return disciplined convex programming result, obtained by calling [CVXR::solve()]
#' 
#' @author Turner Silverthorne
#' @export
run_cvxr_power=function(prob,opts){
  #start    = Sys.time()
  #TimeLimit  PreSolve=0,
  #TODO: configure PreSolve options instead of always using default
  result   = CVXR::solve(prob,verbose=opts$verbose,num_iter=1e9,
                         TimeLimit=opts$time_limit,MIPGapAbs=opts$MIPGapAbs)
  #end      = Sys.time()
  return(result)
}

