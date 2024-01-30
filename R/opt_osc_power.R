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
opt_osc_power=function(dvar0,freqs,Amp,control,
                       nlattice_opts=NULL,alpha=.05,tau=NULL){
if (control$costfun_choice=='svdpower'){
  xout=solve_svdpower(dvar0,freqs,Amp,control,alpha)
  return(xout)
}else if(control$costfun_choice=='svdpower_2lattice'){
  xout=solve_svdpower_2lattice(dvar0,freqs,Amp,control,alpha)
  return(xout)
}else if(control$costfun_choice=='svdpower_discrete'){
  xinds = dvar0
  xout  = solve_svdpower_discrete(xinds,tau,freqs,Amp,control,alpha)
  return(xout)
}else if(control$costfun_choice=='svdpower_2lattice__discrete'){
  if(control$optim_method=='simul_anneal'){
  }else{stop('unknown control$optim_method')}
}else if(control$costfun_choice=='nlattice_power_discrete'){
  if(control$optim_method=='simul_anneal'){
  }else{stop('unknown control$optim_method')}
}else{stop('unknown control$costfun_choice')}
}