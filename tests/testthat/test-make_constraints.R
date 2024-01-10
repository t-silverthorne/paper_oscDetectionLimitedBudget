library(CVXR)
test_that("unrecognized lattice cstr error", {
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'cfun'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  opts$lattice_cstr = 'foo'
  expect_error(make_constraints(x,Lmat,bvec,opts),'opts')
})


test_that('no lattice, length + class of constraints matches ref',{
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'none'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  cstr_ref = list(sum(x)==Constant(opts$Nmeas))
  cstr = make_constraints(x,Lmat,bvec,opts)
  expect_equal(length(cstr),length(cstr_ref))
  expect_equal(class(cstr),class(cstr_ref))
})

test_that('with cfun lattice, length + class of constraints matches ref',{
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'cfun'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  cstr = make_constraints(x,Lmat,bvec,opts)
 
  Nm  = Constant(opts$Nmeas) 
  w   = Constant(bvec)
  mla = Constant(opts$max_lat_active) # convert to CVXR class
  cstr_ref = list( sum(x) <= mla, t(w)%*%x == Nm) #TODO: decide if ineq or eq is easier 
  expect_equal(length(cstr),length(cstr_ref))
  expect_equal(class(cstr),class(cstr_ref))
})

test_that('with lineq lattice, length + class of constraints matches ref',{
  opts              = make_default_opts(solver_type = 'cvxr')
  opts$lattice_cstr = 'lineq'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  cstr = make_constraints(x,Lmat,bvec,opts)
  
  Nm  = Constant(opts$Nmeas) 
  L   = Constant(Lmat) # convert to CVXR class
  b   = bvec-1         # shift trick necessary for constraint
  b   = Constant(b)    # convert to CVXR class
  mla = Constant(opts$max_lat_active) # convert to CVXR class
  cstr_ref = list(sum(x)==Nm,sum(pos(L%*%x - b))<=mla) #TODO: decide if ineq or eq is easier
  
  expect_equal(length(cstr),length(cstr_ref))
  expect_equal(class(cstr),class(cstr_ref))
})

#TODO: check computation manually for easy example

