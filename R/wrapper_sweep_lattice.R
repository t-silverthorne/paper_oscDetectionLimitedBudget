wrapper_sweep_lattice=function(optim_routine,Nvals,Nmeas,cts_flag,dvar0,...){
  fval_best = Inf
  res_best  = NULL
  N1_best   = NULL
  time_tot  = 0
  all_fvals = NULL 
  for (n1 in Nvals){
    n2=Nmeas-n1
    if (cts_flag){
      dvar0[['lat1']]= c(1:n1)/n1 - 1/n1
      dvar0[['lat2']]= c(1:n2)/n2 - 1/n2
      res_loc = optim_routine(dvar0,...)
    }else{
      dvar0$N1=n1
      dvar0$N2=n2
      res_loc = optim_routine(dvar0,...)
    }
    time_tot=time_tot+as.numeric(res_loc$runtime)
    all_fvals=rbind(all_fvals,data.frame(fval=res_loc$fvalue,N1=n1))
    if (res_loc$fvalue<fval_best){
      res_best  = res_loc
      N1_best   = n1
      fval_best = res_loc$fvalue
    }
  }
  return(list(res_best=res_best,N1_best=N1_best,
              time_tot=time_tot,
              all_fvals=all_fvals))
}