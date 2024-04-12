require(gurobi)
require(dplyr)

Nmeas_vals=c(16:48)
w_reg_vals=c(1e-5,1)
drts_vals =c(Inf,1*6,3*6)

setgs =expand.grid(Nmeas=Nmeas_vals,w_reg=w_reg_vals,drts=drts_vals)

# discretization parameters
st_idx       = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
Nfine        = 144 
Nfreq        = 47 
Nmeas        = setgs[st_idx,]$Nmeas #as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
fmin         = 1 
fmax         = 24  # 24
tau          = c(1:Nfine)/Nfine -1/Nfine 
fvec         = seq(from=fmin,to=fmax,length.out=Nfreq)
threads_glob = Sys.getenv("SLURM_CPUS_PER_TASK") 
tlim_glob    = 60*60*2
out_loc_glob = 'figs_for_paper/gurobi_reg_spt/' 
drts         = setgs[st_idx,]$drts
w_reg        = setgs[st_idx,]$w_reg

model=list()

lin_cstr_mode = 'single' # multi or single
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
  # support constraints
  if (lin_cstr_mode=='multi'){
    for (ii in c(1:length(rt_inds))){
      qc_loc              = zero_vec
      qc_loc[rt_inds[ii]] = 1
      lin_csts[[ii+1]]    = qc_loc
      rhs_list[[ii+1]]    = 1
      sense_list[[ii+1]]  = '>'
    }
  }else if(lin_cstr_mode=='single'){
    lin_cstr_vec          = rep(0,Nfine+1)
    rt_inds               = seq(1,Nfine,drts)
    lin_cstr_vec[rt_inds] = 1
    lin_csts[[length(lin_csts)+1]] = lin_cstr_vec  
    rhs_list[[2]]              = sum(lin_cstr_vec)
    sense_list[[2]]            = '>'
  }else{
    stop('unknown lin_cstr_mode')
  }
}

# construct gurobi model
model$A          = lin_csts %>% unlist() %>% matrix(byrow=T,ncol=Nfine+1)
model$sense      = sense_list %>% unlist() 
model$rhs        = rhs_list %>% unlist() 
model$obj        = c(rep(0,Nfine),1)
model$modelsense = 'min'
model$modelname  = 'convex_cosinor_power'
model$vtype      = c(rep('B',Nfine),'C')

# TODO: decide if need to normalize by prior support
# construct matrix for regularization of ncp 
dtau_mat          = outer(tau,tau,'-')
tau_prod_mat      = outer(tau,tau,'*')
ncp_reg_mat       = ((sin(4*pi*fmax*dtau_mat) - sin(4*pi*fmin*dtau_mat))*(4*tau_prod_mat*pi^2 + 1))/(4*dtau_mat)
diag(ncp_reg_mat) = pi*(4*pi^2*tau^2 + 1)*(fmax - fmin)

# quadratic constraints
model$quadcon=list()
for (ii in c(1:Nfreq)){
  freq  = fvec[ii]
  cvec  = matrix(cos(2*pi*freq*tau),nrow=Nfine)
  svec  = matrix(sin(2*pi*freq*tau),nrow=Nfine)
  
  a11   = cvec*cvec
  a22   = svec*svec
  a12   = cvec*svec
  
  qcmat  = a11%*%t(a11) + a22%*%t(a22) +4*a12%*%t(a12)-a11%*%t(a22) - a22%*%t(a11) + w_reg*ncp_reg_mat
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
                        'wreg_',w_reg,
                        '_drts_',drts,
                        '_Nmeas_',Nmeas,'.RDS'))