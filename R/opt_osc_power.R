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
#'  same as previous, in additionlpha
#'  control options necessary for specifying transition function 
#' @author Turner Silverthorne
opt_osc_power=function(dvar0=NULL,freqs,control,
                       nlattice_opts=NULL,tau=NULL,...){
require(lubridate)
require(stringr)
start_time=Sys.time()
if (control$costfun_choice=='svdpower'){
  xout       = solve_svdpower(dvar0,freqs=freqs,control=control,...)
  fvalue     = xout$value
  mtvalue    = xout$par
  xindsvalue = NULL
}else if(control$costfun_choice=='svdpower_2lattice'){
  xout       = solve_svdpower_2lattice(dvar0,freqs=freqs,control=control,...)
  fvalue     = xout$value
  mtvalue    = convert_2lattice_to_state(shift1=dvar0$x0[['shift1']],shift2=xout$par[1],
                                         scale1=xout$par[2],scale2=xout$par[3],dvar0$lat1,dvar0$lat2) 
  xindsvalue = xout$par 
}else if(control$costfun_choice=='svdpower_discrete'){
  xinds      = dvar0
  xout       = solve_svdpower_discrete(xinds=xinds,tau=tau,freqs=freqs,
                                       control=control,...)
  if (control$optim_method=='cvxr'){
    mtvalue    = tau[as.logical(xout[[1]]>1-1e-6)] # TODO: better way of catching this 
    fvalue     = -costfun_svdpower(mt=mtvalue,freqs=freqs,...) 
    xindsvalue = xout[[1]]
  }else if(control$optim_method=='simul_anneal'){
    fvalue     = xout$value
    mtvalue    = tau[xout$par]
    xindsvalue = as.numeric(c(1:control$Nfine) %in% xout$par)
  }else{
    stop('unknown control$optim_method')
  }
}else if(control$costfun_choice=='svdpower_2lattice_discrete'){
  xout=solve_2lattice_svdpower_discrete(dvar0=dvar0,freqs=freqs,
                                        tau=tau,control=control,...)
  fvalue     = xout$value 
  mtvalue    = tau[xinds_from_lat1lat2_pars(N1=dvar0$N1,N2=dvar0$N2,
                                        dx1=xout$par[1],dx2=xout$par[2],
                                        xshift2=xout$par[3])]
  xindsvalue = xout$par 
}else if(control$costfun_choice=='auglattice'){
  xout       = solve_auglattice(dvar0=dvar0,freqs=freqs,control=control,...) 
  fvalue     = xout$value
  mtvalue    = helper_auglattice_to_state(N1=xout$par[1],N2=xout$par[2],
                                     shift2=xout$par[3],scale2=xout$par[4])
  xindsvalue = xout$par
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