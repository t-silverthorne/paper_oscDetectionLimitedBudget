library(CVXR)
test_that("same size var for lineq or none", {
  opts1=make_default_opts(solver_type='cvxr')
  opts2=make_default_opts(solver_type='cvxr')
  opts1$lattice_cstr='none'
  opts2$lattice_cstr='lineq'
  expect_equal(dim(make_variable(opts1)),dim(make_variable(opts2)))
})

test_that('correct class',{
  opts=make_default_opts(solver_type = 'cvxr')
  cc=class(make_variable(opts))
  expect_equal(cc[[1]],'Variable')
})

test_that('unrecgonized lattice options',{
  opts=list(lattice_cstr='foo')
  expect_error(make_variable(opts),'Unrec')
})

test_that('costfun variable is correct size',{
  opts=make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr='cfun'
  Lmat=make_ineqmat(opts)
  var=make_variable(opts,Lmat)
  expect_equal(max(dim(var)),dim(Lmat)[2])
})