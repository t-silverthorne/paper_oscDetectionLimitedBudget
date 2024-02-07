#test_that("function evaluation", {
#  freqs=seq(from=1,to=24,length.out=24)
#  control=list(costfun_choice  = 'nlattice_power_discrete',
#               lattice_cstr    = 'sa_lattice',
#               optim_method    = 'simul_anneal',
#               tfun_choice     = 'uniform',
#               enforce_overlap = 'use-reject',
#               Nmeas           = 16,
#               max_lat = 12, min_lat=5, min_active_lats=2,max_active_lats=4,
#               trace           = 6,
#               REPORT          = 1,
#               maxit           = 50,
#               Nfine=2^7)
#  sa_propfunction(control)
#  #opt_osc_power(freqs=freqs,control=control)
#  
#})
