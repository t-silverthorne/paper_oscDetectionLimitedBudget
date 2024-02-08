solve_svdpower_2lattice=function(dvar0,freqs,Amp,control,alpha,...){
  #unpack
  x0     = dvar0[['x0']]
  lat1   = dvar0[['lat1']]
  lat2   = dvar0[['lat2']]
   
  shift1 = x0[['shift1']]
  shift2 = x0[['shift2']]
  scale1 = x0[['scale1']]
  scale2 = x0[['scale2']]
  
  x0 = c(shift2,scale1,scale2) 
  cfun=function(x){1-costfun_2lattice_svdpower(shift1=0,shift2=x[1],scale1=x[2],
                                               scale2=x[3],lat1,lat2,freqs,Amp,
                                               alpha=alpha,...)}
  if(control$optim_method=='L-BFGS-B'){
    xopt=stats::optim(x0,fn=cfun,gr=NULL,
                 lower=rep(0,length(x0)),
                 upper=rep(1,length(x0)),
                 method='L-BFGS-B',
                 control=list(trace  = control$trace,
                              REPORT = control$REPORT,
                              maxit  = control$maxit))
    return(xopt)
  }else if(control$optim_method=='simul_anneal'){
    tfun=function(x){
      raw=tfun_2lattice_svdpower(shift1=0,shift2=x[1],scale1=x[2],scale2=x[3],
                             lat1,lat2,control)
      return(c(raw$shift2,raw$scale1,raw$scale2))}
    xopt=stats::optim(x0,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
    return(xopt)
  }else{stop('unknown control$optim_method')}
}
