test_that("initializes partition with correct number of points", {
  psize=sample(c('small','medium','large'),1)
  opts=make_default_opts(prob_size=psize,solver_type = 'simulanneal')
  opts$lattice_cstr='none'
  g=sa_propfunction(opts,NULL) 
  nn=g%>% unlist() %>% sum()
  expect_equal(nn,opts$Nmeas)
})

#TODO: check init for lattice constraints