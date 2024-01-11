library(CVXR)
test_that('error handling for run_sa_power',{
  opts = make_default_opts(prob_size='small',
                           solver_type='cvxr')
  opts$costfun_type = 'L1'
  Amat = make_quadmats(opts)
  expect_error(run_sa_power(Amat,opts),
               'You are running simulated annealing, Aquad')
  
  opts$costfun_type = 'Linfty'
  Amat = make_quadmats(opts)
  expect_error(run_sa_power(Amat,opts),
               'You are running simulated annealing, first entry of')
}
)

test_that("sa_power fun eval L1 and Linfty", {
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  
  cft = c('L1','Linfty')
  
  for (cft in cft){
    opts$num_iter = 1e1 
    opts$enforce_overlap = 'ignore' 
    opts$costfun_type = cft 
    Amat = make_quadmats(opts)
    res1 = run_sa_power(Amat,opts)
    res1 = res1$cfunval
    
    opts$lattice_cstr    = 'sa_lattice'
    opts$max_active_lats = 4
    opts$min_active_lats = 1
    opts$min_lat         = 1
    opts$max_lat         = 8 
    res2 =run_sa_power(Amat,opts)
    res2=res2$cfunval 
    expect_gte(res1,0)
    expect_gte(res2,0)
  }
}
)

#TODO: Too expensive to include for now
#test_that("lattice cstr < no cstr", {
#  opts = make_default_opts(prob_size='small',
#                           solver_type='simulanneal')
#  opts$num_iter = 1e4 
#  opts$enforce_overlap = 'ignore' 
#  opts$costfun_type = 'L1'
#  Amat = make_quadmats(opts)
#  res1 = run_sa_power(Amat,opts)
#  res1 = res1$cfunval
#  
#  opts$lattice_cstr    = 'sa_lattice'
#  opts$max_active_lats = 4
#  opts$min_active_lats = 1
#  opts$min_lat         = 1
#  opts$max_lat         = 8 
#  res2 =run_sa_power(Amat,opts)
#  res2=res2$cfunval 
#  expect_lte(res1,res2)
#}
#)

