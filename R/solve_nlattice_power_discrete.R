solve_nlattice_power_discrete=function(Nfine,freqs,Amp=1,control,alpha=.05,...){
  if(control$optim_method=='simul_anneal'){
    cfun=function(xlist){-costfun_nlattice_power_discrete(xlist,Nfine,freqs,Amp,alpha,...)}
    tfun=function(xlist){tfun_nlattice_power_discrete(control,xlist)}
    x0=sa_propfunction(control)
    xout=stats::optim(x0,fn=cfun,gr=tfun,
                      method='SANN',
                      control=list(trace  = control$trace,
                                   REPORT = control$REPORT,
                                   maxit  = control$maxit))
    return(xout)
  }else{stop('uknown control$optim_method')}
}