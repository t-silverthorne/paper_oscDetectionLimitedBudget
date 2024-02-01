solve_svdpower_discrete=function(xinds,tau,freqs,Amp,control,alpha){
  cfun=function(xinds){1-costfun_svdpower_discrete(xinds,tau,freqs,Amp,alpha=alpha)}
  cfun(xinds)
  if(control$optim_method=='simul_anneal'){
    tfun=function(xinds){tfun_svdpower_discrete(xinds,length(tau),control)}
    xopt=stats::optim(xinds,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
    return(xopt)
  }else if(control$optim_method=='cvxr'){
   #TODO 
  }else{stop('unknown control$optim_method')}
}