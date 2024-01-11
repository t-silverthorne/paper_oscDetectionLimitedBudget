library(CVXR)
test_that("lattice addition works", {
  x1  = c(1,0,1,0,0)
  x2  = c(0,1,1,0,1)
  xref= c(1,1,2,0,1)
  x=list(x1,x2)
  opts = list(lattice_cstr='sa_lattice',Nfine=length(x1))
  expect_equal(sa_lat2state(x,opts),xref)
})

test_that('require sa_lattice for opts$lattice_cstr',{
  opts=list(lattice_cstr='foo')
  expect_error(sa_lat2state(list(),opts),'sa_lat2state should')
})

test_that('require x to be a list',{
  opts=list(lattice_cstr='sa_lattice',Nfine=20)
  x=sample(c(0,1),20,replace=T)
  expect_error(sa_lat2state(x,opts),'Input lattice should')
})
  #opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  #opts$lattice_cstr='sa_lattice'
  #opts$max_active_lats = 5
  #opts$max_lat = opts$Nmeas 
  #x=sa_propfunction(opts,NULL) 
