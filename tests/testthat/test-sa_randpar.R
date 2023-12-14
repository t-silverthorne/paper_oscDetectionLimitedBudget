test_that("check uniform sampling for small n", {
  # expect at most 3% deviation from unif distribution freqs
  n=sample(5:8,1)
  N=1e3
  ctab = replicate(N,{sa_randpar(n)}) %>% 
    lapply(function(x){toString(sort(x))}) %>% 
    unlist() %>% 
    table()
  
  perc_error_from_unif = ctab %>% {100*abs(.-N/length(ctab))/N} 
  
  expect_lt(max(perc_error_from_unif),3)
})

test_that("partitions sum correctly", {
  # all partitions should sum to integer
  nrep = 3
  er=c()
  for (ii in c(1:nrep)){
    n=sample(5:20,1)
    N=5e2
    ctab = replicate(N,{sa_randpar(n)}) %>% lapply(function(x){sum(x)}) 
    er[ii]=sum(abs(ctab %>% unlist() %>% unique()-n)) # check if error
  }
  expect_equal(sum(er),0) # check all errors were 0
})

test_that('error for including parent size without options',{
  opts=NULL
  expect_error(sa_randpar(10,opts,100),'Provided size of parent')
 
  opts=list() 
  opts$lattice_cstr='none'
  expect_error(sa_randpar(10,opts,100),'Provided size of parent')
})

test_that('error for including unrecognized lattice constraint',{
opts=list(min_active_lats = 1,
            max_active_lats = 10, 
            min_lat         = 5,
            max_lat         = 10,
            lattice_cstr    = 'wildcard')

expect_error(sa_randpar(10,opts),'Unrecognized lattice constraint')
})

test_that('warning for redundant options',{
  opts=list(lattice_cstr='none')
  expect_warning(sa_randpar(10,opts),'You provided options to sa_randpar')
})

######

test_that('errors for infeasible lattice constraint',{
  opts=list(min_active_lats = 1,
              max_active_lats = 10, 
              min_lat         = 5,
              max_lat         = 10)
  
  # check error when opts provided and enforced
  opts$lattice_cstr = 'sa_lattice'
  expect_error(sa_randpar(200,opts),'Problem is infeasible: increase')
  expect_error(sa_randpar(1,opts),'Number of points to be partitioned')
})


test_that('errors for infeasible lattice constraint with parent size',{
  opts=list(min_active_lats = 1,
              max_active_lats = 10, 
              min_lat         = 5,
              max_lat         = 10,
              lattice_cstr    = 'sa_lattice')
  parsize = 8 
  expect_lte(length(sa_randpar(10,opts,parsize))+parsize,opts$max_active_lats)
  expect_lte(length(sa_randpar(20,opts)),opts$max_active_lats)
  
  parsize=10
  expect_error(sa_randpar(n=20,opts=opts,parent_size=parsize),
               'No feasible partition: you must increase')
  
  expect_error(sa_randpar(n=3,opts=opts,parent_size=10),
               'Number of points to be partitioned is less than')
  parsize=1
  opts$min_active_lats=3
  opts$min_lat =1
  expect_error(sa_randpar(n=1,opts=opts,parent_size=parsize),
               'No feasible partition: you must decrease')
})

test_that('generated partions obey lattice constraints',{
  opts=list(min_active_lats = 1,
            max_active_lats = 10, 
            min_lat         = 5,
            max_lat         = 10,
            lattice_cstr    = 'sa_lattice')
  expect_equal(sa_randpar(opts$min_lat,opts),opts$min_lat)
})
