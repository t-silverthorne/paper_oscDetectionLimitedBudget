#'Solve power optimization using simulated annealing and a pinned 2-lattice
#'
#' @param x0$N1 number of points in the roots of unity grid
#' @param x0$N2 number of points in secondary grid
#' @param x0$shift2 shift of secondary grid relative to start of experiment
#' @param x0$scale2 scale of secondary grid
#' @param freqs list of frequencies for evaluating ojective function
#' @param control$N1min min number of points allowed in N1 lattice, see [tfun_auglattice]
#' @param control$N1max max number of points allowed in N1 lattice, see [tfun_auglattice]
#' @param control$maxit how many iterations of simulated annealing to run. 
#' @param control$trace determines how verbose simulated annealing should be
#' @param control$REPORT determines how verbose simulated annealing should be
#' @param cfuntype set to \code{'ncp'} if non-centrality parameter should be used
#' or \code{'power'} if power should be used in defining cost function.
#' 
#' @note additional inputs \code{...} are passed to [costfun_auglattice] 
#' 
#' @seealso [stats::optim] in the \code{method='SANN'} section for more information 
#' on how iterations are counted. There may be multiple function calls within each iteration.
#'
#'@return result of running simulated annealing using [stats::optim]
#'
#'@author Turner Silverthorne
#'@export
solve_pin2lat=function(x0,freqs,control,...){
  start_time=Sys.time()
  
  N1_init     = x0[['N1']]
  N2_init     = x0[['N2']]
  shift2_init = x0[['shift2']]
  scale2_init = x0[['scale2']]
  
  x0 = c(N1_init,N2_init,shift2_init,scale2_init)
  cfun = function(x){-costfun_auglattice(N1     = x[1],
                                         N2     = x[2],
                                         shift2 = x[3],
                                         scale2 = x[4],
                                         freqs  = freqs,
                                         ...)}
  tfun = function(x){tfun_auglattice(N1     = x[1],
                                     N2     = x[2],
                                     shift2 = x[3],
                                     scale2 = x[4],
                                     freqs  = freqs,
                                     control=control)}
  xout=stats::optim(x0,fn=cfun,gr=tfun,
                    method='SANN',
                    control=list(trace  = control$trace,
                                 REPORT = control$REPORT,
                                 maxit  = control$maxit))
  end_time = Sys.time()
  fvalue     = -xout$value
  mtvalue    = helper_auglattice_to_state(N1=xout$par[1],N2=xout$par[2],
                                     shift2=xout$par[3],scale2=xout$par[4])
  xindsvalue = xout$par
  
  tstamp   = lubridate::now() %>% toString() %>% str_replace(' ','___')
  
  res_full = list(fvalue     = fvalue,
       mtvalue    = mtvalue,
       xindsvalue = xindsvalue,
       timestamp  = tstamp,
       runtime    = end_time-start_time,
       optim_raw  = xout)
}