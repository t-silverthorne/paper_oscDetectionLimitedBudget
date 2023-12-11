#' Generate default options for convex power optimization
#' 
#' @description
#' \code{make_default_opts()} generates optimization options for small medium or large problems.
#' When specifying to your own problem, some options should be kept 
#' 
#' @param prob_size size of problem, either small, medium, or large
#' @return \code{opts} a list containing the following optimization options:
#' @return \code{solver_type} can be cvxr, simulanneal, (not implemented) ctslocal
#' 
#' related to optimization solver:
#' * \code{opts$solver_type} record user-specified \code{solver_type} for downstream functions
#' * \code{opts$verbose} turn on/off intermediate output from Gurobi (\code{TRUE} by default) 
#' * \code{opts$lattice_cstr} decide if you want to constrain design space to only integer
#'                        combinations of certain sub-lattices  (\code{'none'} by default)
#' * \code{opts$costfun_type} either \code{'L1'} which corresponds to averaging over frequencies
#'                        or \code{'Linfty'} which corresponds to maximizing power at worst case frequency
#' * \code{opts$num_iter} iteration cutoff for GUROBI, value depends on \code{prob_size} input
#' 
#' related to time and frequency discretization
#' * \code{opts$Nfine} number of candidate points to include when discretizing time
#' * \code{opts$Nmeas} size of the measurement budget, acts as a constraint
#' * \code{opts$Nfreq} number of frequencies to use in frequency discretization
#' * \code{opts$fmin}  minimum frequency included in prior (for now only uniform prior allowed)
#' * \code{opts$fmax} maximum frequency included in prior (for now only uniform prior allowed)
#'
#' related to lattice construction (only relevant if \code{use_lattice==T}:
#' \itemize{
#'  \item \code{opts$min_dx} smallest lattice spacing 
#'  \item \code{opts$max_lat_active} maximum number of lattices that can be used in constructing optimal solution
#'  \item \code{opts$verbose maximum} number of lattices that can be used in constructing optimal solution
#' }
#' 
#' @author Turner Silverthorne

make_default_opts = function(prob_size='small',solver_type='simulanneal'){
  opts = list( 
    min_dx         = 1, 
    max_lat_active = 5, 
    verbose        = T, 
    fmin           = 1, 
    fmax           = 24,
    lattice_cstr   = 'none', # none, cfun, lineq, sa_lattice
    costfun_type   = 'L1',
    solver_type    = solver_type
    ) 
  
  if (prob_size=='small'){
    opts$Nfine   = 32 
    opts$Nfreq   = 8 
    opts$Nmeas   = 16 
    opts$min_lat = 4 
    opts$max_lat = 4 
    opts$num_iter= 1e6
  }else if (prob_size=='medium'){
    opts$Nfine   = 144 
    opts$Nfreq   = 16 
    opts$Nmeas   = 30 
    opts$min_lat = 6 
    opts$max_lat = 6 
    opts$num_iter= 1e7
  }else if (prob_size=='large'){
    opts$Nfine   = 288 
    opts$Nfreq   = 32 
    opts$Nmeas   = 30 
    opts$min_lat = 6 
    opts$max_lat = 6 
    opts$num_iter= 1e8
  } else if (prob_size=='partial_test'){
    opts=opts
  }else {
    stop("Unknown problem name, use one of the known names")
  }
 
  #TODO document how min and max number of lattices should be handled differently in simulated annealing 
  if (opts$solver_type=='simulanneal'){
    opts$min_active_lats = 1
    opts$max_active_lats = 'adapt' # not hard coded in case Nmeas changes 
    opts$min_lat         = 4
    opts$max_lat         = 'Nmeas' # not hard coded in case Nmeas changes
  } else if (opts$solver_type=='ctslocal'){
    stop('Requested a solver that has not been implemented yet')
  } else{
    stop('Unknown solver requested')
  }
  return(opts)
}