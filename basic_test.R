require(devtools)
library(gurobi)
library(CVXR)
library(dplyr)
devtools::load_all('.')

opts = make_default_opts(prob_size = 'small')
opts$lattice_cstr='lineq' # lineq cfun none
opts$costfun_type='Linfty' # L1 Linfty
Lmat = make_ineqmat(opts)
x    = make_variable(opts,Lmat)
bvec = make_ineqrhs(Lmat,opts)

csts = make_constraints(x,Lmat,bvec,opts)

Aquad = make_quadmats(Lmat,opts)
prob  = make_problem(x,Aquad,csts,opts)
sol   = solve(prob,verbose=opts$verbose,num_iter=opts$num_iter,MIPGapAbs=1e-2,TimeLimit=10)


# loop for benchmarking
