gc()
library(devtools)
library(gurobi)
library(CVXR)
load_all('.')

test=T
# how many batches of simulated annealing 
if (test){
  ens_size = 10 
} else {
  ens_size = 1e2
}

opts=make_default_opts(prob_size='medium',solver_type='simulanneal')

opts$Nfine        = 12*24 # 5 minute intervals 
opts$Nmeas        = 16
opts$fmin         = 1
opts$fmax         = 24
opts$Nfreq        = 2^10 
opts$costfun_type = 'Linfty'
opts$verbose      = F
if (test){
  opts$num_iter     = 15
}else{
  opts$num_iter     = 1e3
}

sa_sols = replicate(ens_size*opts$Nmeas, {NaN})
sa_sols = matrix(sa_sols,nrow=ens_size,ncol=opts$Nmeas)

Aquad=make_quadmats(opts)

for (ii in c(1:ens_size)){
  print(toString(ii))
  xopt         = run_sa_power(Aquad,opts)
  tau          = c(0:opts$Nfine)/opts$Nfine       
  tau          = tau[1:(length(tau)-1)] 
  mt_opt_sa    = tau[xopt$xval==1]
  sa_sols[ii,] = mt_opt_sa
}
saveRDS(sa_sols,paste0('results/sa_sols_','_Nm_',toString(opts$Nmeas),'.RDS'))