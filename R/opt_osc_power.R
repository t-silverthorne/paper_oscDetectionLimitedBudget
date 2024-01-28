#' (WIP) Maximize worst-case power for cosinor-based rhythm detection
#'  Template function for wrapping all optimization methods in suite.
#' 
#' @author Turner Silverthorne
opt_osc_power(control,nlattice_opts=NULL){
if (control$costfun_choice=='svdpower'){
  if(control$optim_method=='L-BFGS-B'){
    
  }else if(control$optim_method=='simul_anneal'){
    
  }else{stop('unknown control$optim_method')}
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