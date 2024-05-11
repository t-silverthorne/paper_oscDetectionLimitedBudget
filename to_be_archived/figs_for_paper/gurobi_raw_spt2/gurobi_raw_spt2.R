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
tlim_glob    = 60*60*2
out_loc_glob = 'figs_for_paper/gurobi_raw_spt2/'
drts            = 2*6

# construct linear constraint matrix

# support constraints
zero_vec        = rep(0,Nfine+1)
rt_inds         = seq(1,Nfine,drts)
lin_csts        = list()
rhs_list        = list()
sense_list      = list()
lin_csts[[1]]   = c(rep(1,Nfine),0)
rhs_list[[1]]   = Nmeas
sense_list[[1]] = '=' 
for (ii in c(1:length(rt_inds))){
  qc_loc              = zero_vec
  qc_loc[rt_inds[ii]] = 1
  lin_csts[[ii+1]]    = qc_loc
  rhs_list[[ii+1]]    = 1
  sense_list[[ii+1]]  = '>'
}
# construct gurobi model
model=list()
model$A          = lin_csts %>% unlist() %>% matrix(byrow=T,ncol=Nfine+1)
model$sense      = sense_list %>% unlist() 
model$rhs        = rhs_list %>% unlist() 
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
saveRDS(sol,file=paste0(out_loc_glob,'sol_gur_',
                        'drts_',drts,
                        '_Nmeas_',Nmeas,'.RDS'))
