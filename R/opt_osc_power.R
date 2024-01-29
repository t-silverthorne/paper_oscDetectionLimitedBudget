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
opt_osc_power=function(mt0,freqs,Amp,control,nlattice_opts=NULL,alpha=.05){
if (control$costfun_choice=='svdpower'){
  xout=solve_svdpower(mt0,freqs,Amp,control,alpha)
}else if(control$costfun_choice=='svdpower_2lattice'){
  if(control$optim_method=='L-BFGS-B'){
      
  }else if(control$optim_method=='simul_anneal'){
    
  }else{stop('unknown control$optim_method')}

}else if(control$costfun_choice=='svdpower_discrete'){

}else if(control$costfun_choice=='svdpower_2lattice__discrete'){
  if(control$optim_method=='simul_anneal'){
  }else{stop('unknown control$optim_method')}
}else if(control$costfun_choice=='nlattice_power_discrete'){
  if(control$optim_method=='simul_anneal'){
  }else{stop('unknown control$optim_method')}
}else{stop('unknown control$costfun_choice')}
}