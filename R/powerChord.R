# methods for Power-CHORD solver
powerChord=function(Nmeas,fmin=1,fmax=24,
                     drts=Inf,w_reg=1,
                     Nfine=144,Nfreq=47,
                     num_threads=1,tlim=60,
                     w_wc=1,MIPGap=1e-5){
   
  model=list()
 
  # construct linear constraints 
  LC = getPcLinCsts(Nfine,drts)
  
  # construct quadratic constraints
  QC = getPcQuadCsts(Nfine,drts)
  
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
  model$quadcon=QC
  
  # gurobi settings
  params = list(TimeLimit=tlim,MIPGap=MIPGap,Presolve=2,
                MIPFocus=3,NumericFocus=3,Threads=num_threads)
  
  # solve
  sol=gurobi(model,params)
  
return(sol)
}