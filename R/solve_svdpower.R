solve_svdpower=function(mt0,freqs,Amp,control,alpha,...){
  cfun=function(mt){-costfun_svdpower(mt,freqs,Amp,alpha=alpha,...)}
  if(control$optim_method=='L-BFGS-B'){
    xopt=stats::optim(mt0,fn=cfun,gr=NULL,
                 lower=rep(0,length(mt0)),
                 upper=rep(1,length(mt0)),
                 method='L-BFGS-B',
                 control=list(trace  = control$trace,
                              REPORT = control$REPORT,
                              maxit  = control$maxit))
    return(xopt)
  }else if(control$optim_method=='simul_anneal'){
    tfun=function(mt){tfun_svdpower(mt,control)}
    xopt=stats::optim(mt0,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
    return(xopt)
  }else{stop('unknown control$optim_method')}
}