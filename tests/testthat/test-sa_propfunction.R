library(CVXR)
test_that("no lattice constraints, initalizes correct number of points", {
  psize=sample(c('small','medium','large'),1)
  opts=make_default_opts(prob_size=psize,solver_type = 'simulanneal')
  opts$lattice_cstr='none'
  g=sa_propfunction(opts,NULL) 
  nn=g%>% unlist() %>% sum()
  expect_equal(nn,opts$Nmeas)
})

test_that('no lattice constraints, transition state has correct number of points',{
  psize=sample(c('small','medium','large'),1)
  opts=make_default_opts(prob_size=psize,solver_type = 'simulanneal')
  opts$lattice_cstr='none'
  x=sa_propfunction(opts,NULL) 
  y=sa_propfunction(opts,x)
  expect_equal(sum(y),opts$Nmeas)
})

test_that('no lattice constraints, transition state differs by only one',{
  psize=sample(c('small','medium','large'),1)
  opts=make_default_opts(prob_size=psize,solver_type = 'simulanneal')
  opts$lattice_cstr='none'
  x=sa_propfunction(opts,NULL) 
  y=sa_propfunction(opts,x)
  expect_equal(length(which(x-y == 1 )),1)
  expect_equal(length(which(x-y == -1 )),1)
  expect_equal(length(which(x-y == 0 )),opts$Nfine-2)
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


test_that('transition states have correct number of points (w lattice)',{
   opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
   opts$lattice_cstr='sa_lattice'
   opts$min_active_lats = 2
   opts$max_active_lats = 5
   opts$min_lat =  4
   opts$max_lat = opts$Nmeas 
   x=sa_propfunction(opts,NULL)
   y=sa_propfunction(opts,x)
   expect_equal(y %>% unlist() %>% sum(),opts$Nmeas)
})

test_that('points get repartitioned correctly',{
  opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = 2
  opts$max_active_lats = 5
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas 
  
  safe_lattice =F
  while (safe_lattice==F){ # diff/intersect only work if all sub lattices are distinct
    x=sa_propfunction(opts,NULL)
    y=sa_propfunction(opts,x)
    safe_lattice = (length(x)==length(unique(x)) & length(y)==length(unique(y)))
  }
 
  xy    = base::intersect(x,y) 
  xminy = base::setdiff(x,y)
  yminx = base::setdiff(y,x)
 
  nshare = xy %>% unlist() %>% sum() 
 
  expect_equal(nshare + (xminy %>% unlist() %>% sum()) ,opts$Nmeas)
  expect_equal(nshare + (yminx %>% unlist() %>% sum()) ,opts$Nmeas)
})


test_that('error for unknown overlap enforcement choice',{
  opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = 2
  opts$max_active_lats = 5
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  x=sa_propfunction(opts,NULL)
  opts$enforce_overlap='wildcard'
  expect_error(sa_propfunction(opts,NULL),'opts')
  expect_error(sa_propfunction(opts,x),'opts')
})


test_that('overlap enforced at initialization',{
  opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = 2
  opts$max_active_lats = 5
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  opts$enforce_overlap='use-reject'
  x=sa_propfunction(opts,NULL)
  no_overlap = Reduce('+',x) %>% {. %in% c(0,1)} %>% prod() %>% as.logical()
  expect_true(no_overlap)
})

test_that('overlap enforced for transition state',{
  opts=make_default_opts(prob_size='medium',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = 2
  opts$max_active_lats = 5
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  opts$enforce_overlap='use-reject'
  x=sa_propfunction(opts,NULL)
  y=sa_propfunction(opts,x)
  no_overlap = Reduce('+',y) %>% {. %in% c(0,1)} %>% prod() %>% as.logical()
  expect_true(no_overlap)
})


test_that('lattice constraints enforced for initial state',{
  # another test does this too, this one is more concise
  opts=make_default_opts(prob_size='large',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = sample(c(1:2),1) 
  opts$max_active_lats = sample(c(opts$min_active_lats+1,opts$min_active_lats+5),1)
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  opts$enforce_overlap='use-reject'
  x=sa_propfunction(opts,NULL)
  expect_lte(length(x),opts$max_active_lats)
  expect_gte(length(x),opts$min_active_lats)
  expect_lte(x %>% lapply(function(t){sum(t)}) %>% unlist() %>% max(), opts$max_lat)
  expect_gte(x %>% lapply(function(t){sum(t)}) %>% unlist() %>% min(), opts$min_lat)
})

test_that('lattice constraints enforced for transition state',{
  opts=make_default_opts(prob_size='large',solver_type = 'simulanneal')
  opts$lattice_cstr='sa_lattice'
  opts$min_active_lats = sample(c(1:2),1) 
  opts$max_active_lats = sample(c(opts$min_active_lats+1,opts$min_active_lats+5),1)
  opts$min_lat =  4
  opts$max_lat = opts$Nmeas
  opts$enforce_overlap='use-reject'
  x=sa_propfunction(opts,NULL)
  y=sa_propfunction(opts,x)
  expect_lte(length(y),opts$max_active_lats)
  expect_gte(length(y),opts$min_active_lats)
  expect_lte(y %>% lapply(function(t){sum(t)}) %>% unlist() %>% max(), opts$max_lat)
  expect_gte(y %>% lapply(function(t){sum(t)}) %>% unlist() %>% min(), opts$min_lat)
})

