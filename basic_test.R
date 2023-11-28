require(devtools)
devtools::load_all('.')

opts = make_default_opts(prob_size = 'small')
opts$lattice_cstr='lineq'
x    = make_variable(opts)
Lmat = make_ineqmat(opts)
bvec = make_ineqrhs(Lmat,opts)

csts = make_constraints(x,Lmat,bvec,opts)

Aquad = make_quadmats(Lmat,opts)
prob  = make_problem(x,Aquad,Lmat,bvec,opts)
sol   = solve_problem(prob,opts)

