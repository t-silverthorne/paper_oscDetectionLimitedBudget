require(devtools)
library(gurobi)
library(CVXR)
library(dplyr)
devtools::load_all()

## CVX programming approach
opts = make_default_opts(prob_size = 'small')
opts$lattice_cstr='none' # lineq cfun none
opts$costfun_type='L1' # L1 Linfty
Lmat = make_ineqmat(opts) # should be depracated
x    = make_variable(opts,Lmat) 
bvec = make_ineqrhs(Lmat,opts)

csts = make_constraints(x,Lmat,bvec,opts)

Aquad = make_quadmats(opts,Lmat)
prob  = make_problem(x,Aquad,csts,opts)

#TODO: Aquad type conversion, wrap this so that you can check types
sol   = solve(prob,verbose=opts$verbose,num_iter=opts$num_iter,MIPGapAbs=1e-2,TimeLimit=10)

devtools::load_all()
opts$num_iter=1e5
run_sa_power(Aquad,opts)

