#test_that("Convex programming warning", {
#  # should warn that ineqmat is not to be used for convex programming
#  opts = make_default_opts(prob_size='small')
#  expect_warning(make_ineqmat(opts),
#                 'The output of this function is only useful for enforcing lattice constraints.*')
#})

test_that('Lmat has more rows than columns',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$min_lat = 3
  opts$max_lat = 4
  opts$lattice_cstr='lineq'
  Lmat=make_ineqmat(opts)
  expect_gte(dim(Lmat)[1],dim(Lmat)[2])
})


test_that('Correct number of pts in lattices',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$Nfine = 12 
  opts$min_lat = 4
  opts$max_lat = 5
  opts$lattice_cstr='lineq'
  Lmat=make_ineqmat(opts)
  lat_sizes = Lmat %>% apply(1,FUN=function(x){sum(x)}) 
  expect_equal(max(lat_sizes),opts$max_lat)
  expect_equal(min(lat_sizes),opts$min_lat)
})

test_that('Lmat entries are 0 or 1',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$min_lat = 3
  opts$max_lat = 4
  opts$lattice_cstr='lineq'
  vals = make_ineqmat(opts) %>% as.vector() %>% unique()
  expect_setequal(vals,c(0,1))
})

test_that('cfun returns transpose of lineq lattice matrix',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$min_lat = 3
  opts$max_lat = 4
  
  opts$lattice_cstr='lineq'
  Lmat1 = make_ineqmat(opts)
  
  opts$lattice_cstr='cfun'
  Lmat2 = make_ineqmat(opts)
  
  expect_equal(Lmat1,t(Lmat2))
})

test_that('unrecognized lattice_cstr',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$lattice_cstr='foo'
  expect_error(make_ineqmat(opts),'Unrecognized lattice')
})

test_that('lattice_cstr=none gives correct response',{
  opts = make_default_opts(prob_size = 'small',solver_type='cvxr')
  opts$lattice_cstr='none'
  expect_null(make_ineqmat(opts))
})