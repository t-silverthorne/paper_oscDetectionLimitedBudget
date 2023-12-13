test_that("uniform sampling for small n", {
  n=3
  N=1e3
  sorted_parts = replicate(N,{sa_randpar(n)}) %>% lapply(function(x){x}) 
  
  padded = lapply(sorted_parts,function(x){
      x=unlist(c(x,replicate(n-length(x),0)))
  })
  summary(padded)
  expect_lt(2 * 2, 4)
})
