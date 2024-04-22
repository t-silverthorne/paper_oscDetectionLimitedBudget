# methods for Power-CHORD solver
powerChord=function(Nmeas,fmin=1,fmax=24,
                     drts=Inf,w_reg=1,
                     Nfine=144,Nfreq=47,
                     num_threads=1,tlim=60){
  
  tau  = c(1:Nfine)/Nfine -1/Nfine 
  fvec = seq(from=fmin,to=fmax,length.out=Nfreq)
  
  model=list()
  
  # construct linear constraint matrix
  zero_vec        = rep(0,Nfine+1)
  rt_inds         = seq(1,Nfine,drts)
  lin_csts        = list()
  rhs_list        = list()
  sense_list      = list()
  lin_csts[[1]]   = c(rep(1,Nfine),0)
  rhs_list[[1]]   = Nmeas
  sense_list[[1]] = '=' 
  
  # if support constraints are active
  if (drts<Inf){
    lin_cstr_vec          = rep(0,Nfine+1)
    rt_inds               = seq(1,Nfine,drts)
    lin_cstr_vec[rt_inds] = 1
    lin_csts[[length(lin_csts)+1]] = lin_cstr_vec  
    rhs_list[[2]]              = sum(lin_cstr_vec)
    sense_list[[2]]            = '>'
  }
  
  # construct gurobi model
  model$A          = lin_csts %>% unlist() %>% matrix(byrow=T,ncol=Nfine+1)
  model$sense      = sense_list %>% unlist() 
  model$rhs        = rhs_list %>% unlist() 
  model$obj        = c(rep(0,Nfine),1)
  model$modelsense = 'min'
  model$modelname  = 'Power_CHORD'
  model$vtype      = c(rep('B',Nfine),'C')

  dtau_mat          = outer(tau,tau,'-')
  tau_prod_mat      = outer(tau,tau,'*')
  ncp_reg_mat       = ((sin(4*pi*fmax*dtau_mat) - sin(4*pi*fmin*dtau_mat))*(4*tau_prod_mat*pi^2 + 1))/(4*dtau_mat)
  diag(ncp_reg_mat) = pi*(4*pi^2*tau^2 + 1)*(fmax - fmin) 
  
  # first set of quadratic constraints
  model$quadcon=list()
  for (ii in c(1:Nfreq)){
    freq  = fvec[ii]
    cvec  = matrix(cos(2*pi*freq*tau),nrow=Nfine)
    svec  = matrix(sin(2*pi*freq*tau),nrow=Nfine)
    
    a11   = cvec*cvec
    a22   = svec*svec
    a12   = cvec*svec
    
    qcmat  = a11%*%t(a11) + a22%*%t(a22) +4*a12%*%t(a12)-a11%*%t(a22) - a22%*%t(a11)+w_reg*ncp_reg_mat
    qcmat  = cbind(rbind(qcmat,rep(0,Nfine)),rep(0,Nfine+1))
    qc     = c(rep(0,Nfine),-1)
    beta   = 0
    sense  = '<'
    model$quadcon[[ii]] = list(Qc=qcmat,q=qc,rhs=beta,sense=sense)
  }
  
   
  # gurobi settings
  params = list(TimeLimit=tlim,MIPGap=1e-5,Presolve=2,
                MIPFocus=3,NumericFocus=3,Threads=num_threads)
  
  # solve
  sol=gurobi(model,params)
  
return(sol)
}