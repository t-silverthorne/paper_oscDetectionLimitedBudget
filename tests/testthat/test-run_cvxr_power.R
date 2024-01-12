library(CVXR)
test_that("cvxr for L1 Linfty runs, no lattice cstr", {
  opts              = make_default_opts(solver_type = 'cvxr',
                                        prob_size = 'small')
  opts$lattice_cstr = 'none'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  csts              = make_constraints(x,Lmat,bvec,opts)
  prob              = make_problem(x,Aquad,csts,opts)
  
  opts$costfun_type = 'L1'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-2
  run_cvxr_power(prob,opts)
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
  
  opts$costfun_type = 'Linfty'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-2
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
})

test_that("cvxr for L1 Linfty runs, cfun lattice_cstr", {
  opts              = make_default_opts(solver_type = 'cvxr',
                                        prob_size = 'small')
  opts$lattice_cstr = 'cfun'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  csts              = make_constraints(x,Lmat,bvec,opts)
  prob              = make_problem(x,Aquad,csts,opts)
  
  opts$costfun_type = 'L1'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-1
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
  
  opts$costfun_type = 'Linfty'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-1
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
})


test_that("cvxr for L1 Linfty runs, lineq lattice_cstr", {
  opts              = make_default_opts(solver_type = 'cvxr',
                                        prob_size = 'small')
  opts$lattice_cstr = 'lineq'
  Lmat              = make_ineqmat(opts)
  bvec              = make_ineqrhs(Lmat,opts)
  Aquad             = make_quadmats(opts,Lmat)
  x                 = make_variable(opts,Lmat)
  csts              = make_constraints(x,Lmat,bvec,opts)
  prob              = make_problem(x,Aquad,csts,opts)
  
  opts$costfun_type = 'L1'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-1
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
  
  opts$costfun_type = 'Linfty'
  opts$verbose      = F
  opts$MIPGapAbs    = 1e-1
  res=run_cvxr_power(prob,opts)
  expect_gte(res$value,0)
})