library(CVXR)
test_that("sa_cfunpwr handles expected errors", {
  x    = NULL
  Amat = NULL
  opts=list(solver_type='foo')
  expect_error(sa_cfunpwr(x,Amat,opts),'Wrong solver')
  
  opts=list(solver_type='simulanneal',lattice_cstr = 'foo')
  expect_error(sa_cfunpwr(x,Amat,opts),
               'Unrecognized lattice constraint')
 
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$lattice_cstr = 'none'
  x                 = sa_propfunction(opts)
  opts$costfun_type = 'foo'
  expect_error(sa_cfunpwr(x,Amat,opts),'Cost function')
  
  opts$lattice_cstr = 'sa_lattice'
  opts$max_active_lats = 4
  opts$min_active_lats = 3
  opts$min_lat         = 2
  opts$max_lat         = 8
  opts$enforce_overlap = 'use-reject' 
  x                 = sa_propfunction(opts)
  opts$costfun_type = 'foo'
  expect_error(sa_cfunpwr(x,Amat,opts),'Cost function')
})

test_that('costfunction non-negative',{
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$lattice_cstr = 'sa_lattice'
  opts$max_active_lats = 4
  opts$min_active_lats = 3
  opts$min_lat         = 2
  opts$max_lat         = 8
  opts$enforce_overlap = 'use-reject' 
  x                 = sa_propfunction(opts)
  opts$costfun_type = 'L1'
  Amat = make_quadmats(opts)
  expect_gte(sa_cfunpwr(x,Amat,opts),0) 
  
  opts$costfun_type = 'Linfty'
  Amat = make_quadmats(opts)
  expect_gte(sa_cfunpwr(x,Amat,opts),0) 
})


test_that('detailed errors for erroneously mixing L1/Linfty settings',{
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$lattice_cstr = 'sa_lattice'
  opts$max_active_lats = 4
  opts$min_active_lats = 3
  opts$min_lat         = 2
  opts$max_lat         = 8
  opts$enforce_overlap = 'use-reject' 
  x                 = sa_propfunction(opts)
  
  opts$costfun_type = 'L1'
  Amat = make_quadmats(opts)
  opts$costfun_type = 'Linfty'
  expect_error(sa_cfunpwr(x,Amat,opts),'For Linfty') 
  
  opts$costfun_type = 'Linfty'
  Amat = make_quadmats(opts)
  opts$costfun_type = 'L1'
  expect_error(sa_cfunpwr(x,Amat,opts),'For L1') 
})
#TODO computation for known reference value