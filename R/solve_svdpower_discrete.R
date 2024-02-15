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
    result = CVXR::solve(prob,verbose=control$cvxr_verbose,num_iter=control$maxit,
                         TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,
                         Presolve=2,MIPFocus=3)
  }else{stop('unknown control$optim_method')}
}