library(CVXR)
test_that("global opts always same", {
  # check problem size does not change global options
  opts_ref=make_default_opts(prob_size = 'partial_test','cvxr')
  for (psize in c('small','medium','large')){
    opts=make_default_opts(prob_size = 'small','cvxr')
    expect_equal(list(opts$min_dx,
                      opts$max_lat_active,
                      opts$verbose,
                      opts$fmin,
                      opts$fmax,
                      opts$lattice_cstr,
                      opts$costfun_type),
                 list(opts_ref$min_dx,
                      opts_ref$max_lat_active,
                      opts_ref$verbose,
                      opts_ref$fmin,
                      opts_ref$fmax,
                      opts_ref$lattice_cstr,
                      opts_ref$costfun_type)
                )
  }
})

test_that('unknown problem name gives error',{
  expect_error(make_default_opts(prob_size='huge'),"Unknown problem name, use one of the known names")
})

test_that('CVXR problem sizes nest as expected',{
  stype = 'cvxr' 
  opts_small = make_default_opts(prob_size = 'small',stype)
  opts_med   = make_default_opts(prob_size = 'medium',stype)
  opts_large = make_default_opts(prob_size = 'large',stype)
  
  expect_lte(opts_small$Nfine,opts_med$Nfine)
  expect_lte(opts_small$Nfreq,opts_med$Nfreq)
  expect_lte(opts_small$Nmeas,opts_med$Nmeas)
  expect_lte(opts_small$min_lat,opts_med$min_lat)
  expect_lte(opts_small$max_lat,opts_med$max_lat)
  expect_lte(opts_small$num_iter,opts_med$num_iter)
  
  expect_lte(opts_med$Nfine,opts_large$Nfine)
  expect_lte(opts_med$Nfreq,opts_large$Nfreq)
  expect_lte(opts_med$Nmeas,opts_large$Nmeas)
  expect_lte(opts_med$min_lat,opts_large$min_lat)
  expect_lte(opts_med$max_lat,opts_large$max_lat)
  expect_lte(opts_med$num_iter,opts_large$num_iter)
})