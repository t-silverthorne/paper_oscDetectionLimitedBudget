solve_cvxr_spt= function(control,Threads,use_spt=F,drts=NULL,pads=NULL){
  Aquad = make_quadmats(control)
  x     = Variable(control$Nfine,boolean=T)
 
  #  ensure support distributed throughout domain
  rt_inds=seq(1,control$Nfine,drts)
  csts  = list(sum(x)==control$Nmeas)
  
  if (use_spt){
    for (ii in c(1:length(rt_inds))){
      cst_spt = sum(x[c(rt_inds[ii]:(rt_inds[ii]+pads))])>=1
      csts = append(csts,cst_spt)
    }
  }

  # convert to CVXR datatype
  for (ii in c(1:length(Aquad))){
    Aquad[[ii]]=Constant(Aquad[[ii]])
  }
  
  str_prefix='prob=Problem(Minimize('
  big_str = paste(lapply(1:length(Aquad) %>% as.list(),
         function(ind){
           paste0('quad_form(x,Aquad[[',ind,']])')
         }),collapse=',')
 
  str_suffix=paste0('max_elemwise(',big_str,')),csts)')
 
  strp=paste0(str_prefix,str_suffix)
  eval(parse(text=strp))

  xout=CVXR::solve(prob,verbose=control$cvxr_verbose,num_iter=control$maxit,solver="GUROBI",
                       TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,MIPGap=control$MIPGap,
                       Presolve=control$PreSolve,MIPFocus=control$MIPFocus,
                       Threads=Threads,NodefileStart=control$NodefileStart)
  tau = c(1:control$Nfine)/control$Nfine - 1/control$Nfine
  mtvalue    = tau[as.logical(xout[[1]]>1-1e-6)] # TODO: better way of catching this 
  
  return(mtvalue)
}