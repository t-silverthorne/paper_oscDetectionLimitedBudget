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