library(CVXR)
test_that("triggered if solver_type wrong", {
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'cfun'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  csts              = make_constraints(x,Lmat,bvec,opts)
  opts$solver_type  = 'foo'
  expect_error(make_problem(x,Aquad,csts,opts),'make_problem')
})

test_that("triggered if costfun_type wrong", {
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'cfun'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  csts              = make_constraints(x,Lmat,bvec,opts)
  opts$costfun_type = 'foo'
  expect_error(make_problem(x,Aquad,csts,opts),'opts')
})

test_that("trigerred if costfun_type not recognized", {
  expect_equal(2 * 2, 4)
})

#TODO: check computation manually