require(devtools)
library(gurobi)
library(CVXR)
library(dplyr)

# load this package 
devtools::load_all('.')

opts = make_default_opts(prob_size = 'small')     # default options
opts$lattice_cstr='cfun'                          # options: lineq cfun none
opts$costfun_type='L1'                        # options:  L1 Linfty

Lmat = make_ineqmat(opts)                         # matrix for lattice constraints
x    = make_variable(opts,Lmat)                   # variable may require Lmat size
bvec = make_ineqrhs(Lmat,opts)                    # rhs of lattice constraints

csts = make_constraints(x,Lmat,bvec,opts)         # list, integer and lattice constraints

Aquad = make_quadmats(Lmat,opts)                  # symmetric matrices for cost fun 

prob  = make_problem(x,Aquad,csts,opts)           # make CVXR problem
sol   = solve(prob,                               # solve with GUROBI
              verbose=opts$verbose,
              num_iter=opts$num_iter,
              MIPGapAbs=1e-2,
              TimeLimit=5)


