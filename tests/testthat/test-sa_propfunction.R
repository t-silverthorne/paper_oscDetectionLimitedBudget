test_that("without lattice constraints, initalizes correct number of points", {
  psize=sample(c('small','medium','large'),1)
  opts=make_default_opts(prob_size=psize,solver_type = 'simulanneal')
  opts$lattice_cstr='none'
  g=sa_propfunction(opts,NULL) 
  nn=g%>% unlist() %>% sum()
  expect_equal(nn,opts$Nmeas)
})


test_that('error for unknown lattice constraints',{
  opts=list(lattice_cstr='wildcard')
  expect_error(sa_propfunction(opts),'choice of lattice constraint')
  expect_error(sa_propfunction(opts,x=c(1,0,0,0,1)),'choice of lattice constraint')
})

test_that('correct number of points, lattice constrained state generation',{
  nrep=10
  cvec = replicate(nrep,{0})
  for (ii in c(1:nrep)){
   opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
   opts$lattice_cstr='sa_lattice'
   opts$max_active_lats = 5
   opts$max_lat = opts$Nmeas 
   x=sa_propfunction(opts,NULL)
   cvec[ii] = x %>% unlist() %>% sum()
  }
  expect_equal(cvec,replicate(nrep,opts$Nmeas))
})

test_that('initialized states satisfy lattice constraints',{
  nrep=10
  cvec1 = replicate(nrep,{TRUE})
  cvec2 = replicate(nrep,{0})
  cvec3 = replicate(nrep,{0})
  for (ii in c(1:nrep)){
   opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
   opts$lattice_cstr='sa_lattice'
   opts$max_active_lats = 5
   opts$max_lat = opts$Nmeas 
   x=sa_propfunction(opts,NULL)
   cvec1[ii] = ((length(x) <= opts$max_active_lats) &
                  (length(x) >= opts$min_active_lats))
   cvec2[ii] = x %>% lapply(function(t){sum(t)}) %>% unlist() %>% max() 
   cvec3[ii] = x %>% lapply(function(t){sum(t)}) %>% unlist() %>% min() 
  }
  expect_equal(cvec1,replicate(nrep,T))
  expect_lte(max(cvec2),opts$max_lat)
  expect_gte(max(cvec3),opts$min_lat)
})


#TODO: transition states have correct number of points (w/wo lattice)
#TODO: transition states satisfy lattice constraints (w constraints imposed)
