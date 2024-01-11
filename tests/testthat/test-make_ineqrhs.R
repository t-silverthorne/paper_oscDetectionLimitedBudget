library(CVXR)
test_that('error for unrecognized',{
  opts=list(lattice_cstr='foo')
  expect_error(make_ineqrhs(NULL,opts),'Unrec')
})

test_that("correct ineqrhs val", {
  opts=make_default_opts(solver='cvxr')
  opts$lattice_cstr='cfun'
  Lmat=make_ineqmat(opts)
  expect_equal(make_ineqrhs(Lmat,opts),colSums(Lmat))
  opts$lattice_cstr='lineq'
  Lmat=make_ineqmat(opts)
  expect_equal(make_ineqrhs(Lmat,opts),rowSums(Lmat))
})

test_that('handling of lattice_cstr=none or unrecognized',{
  opts=make_default_opts(solver='cvxr')
  opts$lattice_cstr='none'
  expect_null(make_ineqrhs(NULL,opts))
  
  opts$lattice_cstr='foo'
  expect_error(make_ineqrhs(NULL,opts))
})
