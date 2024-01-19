gc()
library(devtools)
library(gurobi)
library(CVXR)
load_all('.')

test              = T

opts              = make_default_opts(prob_size='medium',solver_type='cvxr')
opts$Nfine        = 12*24 # 5 minute intervals 
opts$Nmeas        = 16
opts$fmin         = 1
opts$fmax         = 24
opts$costfun_type = 'Linfty'
opts$Nfreq        = 2^9 #want to max this before out of memory in pre-solve
opts$verbose      = T
if (test){
  opts$time_limit  = 10 # units of seconds
}else{
  opts$time_limit  = 60*60*24*3 #  3 days 
}

x     = make_variable(opts)
csts  = make_constraints(x,NULL,NULL,opts)

Aquad = make_quadmats(opts)
prob  = make_problem(x,Aquad,csts,opts)
start=Sys.time()
xopt = CVXR::solve(prob,verbose=opts$verbose,num_iter=1e9,
                   TimeLimit=opts$time_limit,
                   MIPGapAbs=opts$MIPGapAbs,
                   MIPFocus=2,
                   Presolve=2
                    )
end=Sys.time()
end-start

tau          = c(0:opts$Nfine)/opts$Nfine       
tau          = tau[1:(length(tau)-1)] 
cvxr_sol     = tau[xopt$getValue(objet = x)==1]

#saveRDS(cvxr_sol,paste0('results/cvxr_sol_','_Nm_',toString(opts$Nmeas),'.RDS'))
#saveRDS(xopt,paste0('results/cvxr_xopt_','_Nm_',toString(opts$Nmeas),'.RDS'))
#saveRDS(opts,paste0('results/cvxr_opts_','_Nm_',toString(opts$Nmeas),'.RDS'))
