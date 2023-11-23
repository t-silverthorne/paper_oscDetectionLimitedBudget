#' Generate default optimization options
#' 
#' @description
#' \code{make_default_opts()} generates optimization options for small medium or large problems.
#' When specifying to your own problem, some options should be kept 
#' 
#' @param prob_size size of problem, either small, medium, or large
#' @return \code{opts} a list containing the following optimization options:
#' 
#' related to optimization solver:
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

#TODO: expand documentation to contain all outputs of opts
make_default_opts = function(prob_size='small'){
  opts = list( 
    min_dx         = 1, 
    max_lat_active = 5, 
    verbose        = T, 
    fmin           = 1, 
    fmax           = 24,
    lattice_cstr   = 'none', # none, cfun, lineq
    costfun_type   = 'L1'
    ) 
  if (prob_size=='small'){
    opts$Nfine   = 64
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
  }else{
    opts = NaN
  }
  return(opts)
}