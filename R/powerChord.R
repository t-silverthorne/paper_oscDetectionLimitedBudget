# methods for Power-CHORD solver
powerChord=function(Nmeas,fmin=1,fmax=24,
                     drts=Inf,w_reg=1,
                     Nfine=144,Nfreq=47,
                     num_threads=1,tlim=60,
                     w_wc=1,MIPGap=1e-5,return_model=F){
  model=list()
 
  # construct linear constraints 
  LC = getPcLinCsts(Nmeas,Nfine,drts)
  
  # construct quadratic constraints
  QC = getPcQuadCsts(fmin,fmax,Nfreq,w_reg,w_wc,Nfine)
  
  # add linear constraints to gurobi model 
  model$A          = LC$A
  model$sense      = LC$sense_list %>% unlist() 
  model$rhs        = LC$rhs_list %>% unlist() 
  model$obj        = c(rep(0,Nfine),1)
  
  # model settings
  model$modelsense = 'min'
  model$modelname  = 'Power_CHORD'
  model$vtype      = c(rep('B',Nfine),'C')

  # add quadratic constraints
  if (w_wc == 0){
    model$quadcon=list()
    model$quadcon[[1]]=QC
  }else{
    model$quadcon=QC
  }
  
  if(return_model){
    return(model) # useful for unit tests
  }else{
    # gurobi settings
    params = list(TimeLimit=tlim,MIPGap=MIPGap,Presolve=2,
                  MIPFocus=3,NumericFocus=3,Threads=num_threads)
    
    # solve
    sol=gurobi(model,params)
    return(sol)    
  }

}