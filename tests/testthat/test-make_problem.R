library(CVXR)
test_that("triggered if wrong costfun", {
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'none'
  Aquad             = make_quadmats(opts)
  x                 = make_variable(opts)
  opts$solver_type  = 'foo'
  expect_error(make_problem(x,Aquad,opts),'unrec')
})

