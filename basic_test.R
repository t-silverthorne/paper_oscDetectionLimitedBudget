require(devtools)
library(gurobi)
library(CVXR)
library(dplyr)
# load this package 
devtools::load_all('.')

## CVX programming approach
opts = make_default_opts(prob_size = 'small')     # small medium large
opts$lattice_cstr='none'                          # lineq cfun none
opts$costfun_type='Linfty'                        # L1 Linfty

Lmat = make_ineqmat(opts)                         # matrix for lattice constraints
x    = make_variable(opts,Lmat)                   # variable may require Lmat size
bvec = make_ineqrhs(Lmat,opts)                    # rhs of lattice constraints

Aquad = make_quadmats(opts,Lmat)
prob  = make_problem(x,Aquad,csts,opts)

#TODO: Aquad type conversion, wrap this so that you can check types
sol   = solve(prob,verbose=opts$verbose,num_iter=opts$num_iter,MIPGapAbs=1e-2,TimeLimit=10)

devtools::load_all()

## simulated annealing approach
opts$num_iter=1e5
run_sa_power(Aquad,opts)

