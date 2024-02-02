test_that("function evaluation", {
  opts=make_default_opts(prob_size='medium',solver_type='simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$max_active_lats=4
  opts$min_active_lats=4
  opts$max_lat=10
  opts$min_lat=4
  control=list(tfun_choice='default')
  control=c(control,opts)
  expect_equal(length(tfun_nlattice_power_discrete(control,NULL)),4)
})
