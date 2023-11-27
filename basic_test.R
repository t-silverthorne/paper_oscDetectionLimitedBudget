require(devtools)
devtools::load_all('.')

opts = make_default_opts(prob_size = 'small')
opts$lattice_cstr='lineq'
x    = make_variable(opts)
Lmat = make_ineqmat(opts)
bvec = make_ineqvec(Lmat,opts)

make_constraints(x,Lmat,bvec,opts)