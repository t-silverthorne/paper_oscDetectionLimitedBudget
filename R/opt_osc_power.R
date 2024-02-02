#' (WIP) Maximize worst-case power for cosinor-based rhythm detection
#'  Template function for wrapping all optimization methods in suite.
#' 
#' @param control$costfun_choice
#' @param control$optim_method
#' 
#' svdpower with L-BFGS-B requires:
#'  control$trace
#'  control$REPORT
#'  control$maxit
#' svdpower with simul_anneal requires:
#'  same as previous, in addition
#'  control options necessary for specifying transition function 
#' @author Turner Silverthorne
opt_osc_power=function(dvar0=NULL,freqs,Amp=1,control,
                       nlattice_opts=NULL,alpha=.05,tau=NULL){
require(lubridate)
require(stringr)
start_time=Sys.time()
if (control$costfun_choice=='svdpower'){
  xout       = solve_svdpower(dvar0,freqs,Amp,control,alpha)
  # output handling
  fvalue     = xout$value
  mtvalue    = xout$par
  xindsvalue = NULL
}else if(control$costfun_choice=='svdpower_2lattice'){
  xout       = solve_svdpower_2lattice(dvar0,freqs,Amp,control,alpha)
  # output handling
  fvalue     = xout$value
  mtvalue    = convert_2lattice_to_state(shift1=dvar0$x0[['shift1']],shift2=xout$par[1],
                                         scale1=xout$par[2],scale2=xout$par[3],dvar0$lat1,dvar0$lat2) 
  xindsvalue = xout$par 
}else if(control$costfun_choice=='svdpower_discrete'){
  xinds      = dvar0
  xout       = solve_svdpower_discrete(xinds,tau,freqs,Amp,control,alpha)
  fvalue     = xout$value
  mtvalue    = NULL#convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2) 
  xindsvalue = NULL#xout$par 
}else if(control$costfun_choice=='svdpower_2lattice_discrete'){
  xout=solve_2lattice_svdpower_discrete(dvar0,freqs,tau,Amp,control,alpha)
  fvalue     = xout$value
  mtvalue    = NULL#convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2) 
  xindsvalue = NULL#xout$par 
}else if(control$costfun_choice=='nlattice_power_discrete'){
  stop("Currently not able to handle n-lattice constraints. No simulated annealing routine implemented")
  #xout=solve_nlattice_power_discrete(Nfine,freqs,Amp,control,alpha)
  #return(xout)
}else{
  stop('unknown control$costfun_choice')
}
  end_time = Sys.time()
  tstamp   = lubridate::now() %>% toString() %>% str_replace(' ','___')
  
  # return optimizer output and optimization settings
  res_full = list(fvalue     = fvalue,
       mtvalue    = mtvalue,
       xindsvalue = xindsvalue,
       timestamp  = tstamp,
       runtime    = end_time-start_time,
       optim_raw  = xout)
  return(res_full)
}