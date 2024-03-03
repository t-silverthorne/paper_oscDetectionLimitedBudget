solve_svdpower=function(mt0,freqs,control,...){
  cfun=function(mt){-costfun_svdpower(mt,freqs,...)}
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

solve_svdpower_discrete=function(xinds,tau,freqs,control,...){
  if(control$optim_method=='simul_anneal'){
    cfun=function(xinds){-costfun_svdpower_discrete(xinds=xinds,tau=tau,freqs=freqs,...)}
    tfun=function(xinds){tfun_svdpower_discrete(xinds,length(tau),control)}
    xopt=stats::optim(xinds,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
    return(xopt)
  }else if(control$optim_method=='cvxr'){
    #if(seq(from=control$fmin,to=control$fmax,length.out=Nfreq)!=freqs){
    #  stop('inconsistent freq specification in cvxr')
    #}
    Aquad             = make_quadmats(control)
    x                 = make_variable(control)
    prob              = make_problem(x,Aquad,control,...)
    result = CVXR::solve(prob,verbose=control$cvxr_verbose,solver="GUROBI",
                         num_iter=control$maxit,
                         TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,
                         Presolve=2,MIPFocus=3)
  }else{stop('unknown control$optim_method')}
}

solve_svdpower_2lattice=function(dvar0,freqs,control,...){
  #unpack
  x0     = dvar0[['x0']]
  lat1   = dvar0[['lat1']]
  lat2   = dvar0[['lat2']]
   
  shift1 = x0[['shift1']]
  shift2 = x0[['shift2']]
  scale1 = x0[['scale1']]
  scale2 = x0[['scale2']]
  
  x0 = c(shift2,scale1,scale2) 
  cfun=function(x){-costfun_2lattice_svdpower(shift1=0,shift2=x[1],scale1=x[2],
                                              scale2=x[3],
                                              lat1=lat1,lat2=lat2,
                                              freqs=freqs,...)}
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


solve_auglattice=function(dvar0,freqs,control,cfuntype='ncp',...){
  N1_init     = dvar0[['N1']]
  N2_init     = dvar0[['N2']]
  shift2_init = dvar0[['shift2']]
  scale2_init = dvar0[['scale2']]
 
  x0 = c(N1_init,N2_init,shift2_init,scale2_init)
  cfun = function(x){-costfun_auglattice(N1     = x[1],
                                         N2     = x[2],
                                         shift2 = x[3],
                                         scale2 = x[4],
                                         freqs  = freqs,
                                         cfuntype=cfuntype,
                                         ...)}
  tfun = function(x){tfun_auglattice(N1     = x[1],
                                     N2     = x[2],
                                     shift2 = x[3],
                                     scale2 = x[4],
                                     freqs  = freqs,
                                     control=control)}
   xopt=stats::optim(x0,fn=cfun,gr=tfun,
                  method='SANN',
                  control=list(trace  = control$trace,
                               REPORT = control$REPORT,
                               maxit  = control$maxit))
   return(xopt)
}

solve_2lattice_svdpower_discrete=function(dvar0,freqs,tau,control,...){
x0      = dvar0[['x0']]
N1      = dvar0[['N1']]
N2      = dvar0[['N2']]
dx1     = x0[['dx1']]
dx2     = x0[['dx2']]
xshift2 = x0[['xshift2']]

if (control$optim_method=='simul_anneal'){
  cfun=function(x){
    if (any(is.na(x))){
      stop('NA value in meastimes')
    }else{
      return(-costfun_2lattice_svdpower_discrete(N1=N1,dx1=x[1],N2=N2,
                                                  dx2=x[2],xshift2=x[3],
                                                  tau=tau,freqs=freqs,...))
    }
  }
  tfun=function(x){
    if (any(is.na(tau[xinds_from_lat1lat2_pars(N1=N1,dx1=x[1],N2=N2,dx2=x[2],xshift2=x[3])]))){
      stop('NA value generated by pars')
    }else{
      raw=tfun_2lattice_svdpower_discrete(N1=N1,dx1=x[1],N2=N2,dx2=x[2],xshift2=x[3],Nfine=length(tau),control=control)
      return(c(raw$dx1,raw$dx2,raw$xshift2))  
    }
    }
  xopt=stats::optim(x0,fn=cfun,gr=tfun,
                method='SANN',
                control=list(trace  = control$trace,
                             REPORT = control$REPORT,
                             maxit  = control$maxit))
  return(xopt)
}else{stop('unknown control$optim_method')}
}
  