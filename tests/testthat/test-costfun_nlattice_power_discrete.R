test_that("function evaluation, alpha handling", {
  opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = 2
  opts$max_active_lats = 5
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  opts$enforce_overlap='use-reject'
  x=sa_propfunction(opts,NULL)
  expect_gte(costfun_nlattice_power_discrete(x,opts$Nfine,c(1,2,3),alpha=.05),
             costfun_nlattice_power_discrete(x,opts$Nfine,c(1,2,3),alpha=.01))
})
