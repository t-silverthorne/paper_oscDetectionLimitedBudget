gc()
library(devtools)
library(gurobi)
library(CVXR)
load_all('.')

test              = F

opts              = make_default_opts(prob_size='medium',solver_type='cvxr')
opts$Nfine        = 12*24 # 5 minute intervals 
opts$Nmeas        = 16
opts$fmin         = 1
opts$fmax         = 24
opts$costfun_type = 'Linfty'
opts$Nfreq        = 2^8 #want to max this before out of memory in pre-solve
opts$verbose      = T
if (test){
  opts$time_limit  = 5 # units of seconds
}else{
  opts$time_limit  = 60*60*24*3 #  3 days 
}

x     = make_variable(opts)
csts  = make_constraints(x,NULL,NULL,opts)
Aquad = make_quadmats(opts)
prob  = make_problem(x,Aquad,csts,opts)
xopt  = run_cvxr_power(prob,opts)

cvxr_sol = tau[xopt$getValue(objet = x)==1]

saveRDS(cvxr_sol,'results/cvxr_sol.RDS')