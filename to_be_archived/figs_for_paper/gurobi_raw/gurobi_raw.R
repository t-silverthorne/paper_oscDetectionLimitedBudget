require(gurobi)
require(dplyr)

# discretization parameters
Nfine = 144      
Nfreq = 47 
Nmeas = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
fmin  = 1
fmax  = 24
tau   = c(1:Nfine)/Nfine -1/Nfine 
fvec  = seq(from=fmin,to=fmax,length.out=Nfreq)
threads_glob = Sys.getenv("SLURM_CPUS_PER_TASK") 
tlim_glob    = 2*60*60
out_loc_glob = 'figs_for_paper/gurobi_raw/'


# construct gurobi model
model=list()
model$A          = matrix(c(rep(1,Nfine),0),nrow=1)
model$sense      = '='
model$rhs        = Nmeas
model$obj        = c(rep(0,Nfine),1)
model$modelsense = 'min'
model$modelname  = 'convex_cosinor_power'
model$vtype      = c(rep('B',Nfine),'C')

# quadratic constraints
model$quadcon=list()
for (ii in c(1:Nfreq)){
  freq  = fvec[ii]
  cvec  = matrix(cos(2*pi*freq*tau),nrow=Nfine)
  svec  = matrix(sin(2*pi*freq*tau),nrow=Nfine)
  
  a11   = cvec*cvec
  a22   = svec*svec
  a12   = cvec*svec
  
  qcmat  = a11%*%t(a11) + a22%*%t(a22) +4*a12%*%t(a12)-a11%*%t(a22) - a22%*%t(a11)
  qcmat  = cbind(rbind(qcmat,rep(0,Nfine)),rep(0,Nfine+1))
  qc     = c(rep(0,Nfine),-1)
  beta   = 0
  sense  = '<'
  model$quadcon[[ii]] = list(Qc=qcmat,q=qc,rhs=beta,sense=sense)
}

# gurobi settings
params = list(TimeLimit=tlim_glob,MIPGap=1e-5,
              Presolve=2,MIPFocus=3,NumericFocus=3,Threads=threads_glob)

# solve
sol=gurobi(model,params)
saveRDS(sol,file=paste0(out_loc_glob,'sol_gur_free_Nmeas_',Nmeas,'.RDS'))
