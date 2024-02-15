test_that("function evaluation", {
  mt = runif(30)
  expect_no_error(deig_dfreq(mt,runif(1)))
})

test_that('agrees with finite difference',{
  set.seed(1)
  mt = runif(20)
  
  freq = 2.5
  df   = 1e-5
  
  dfreq_FD    = (getMinEig(getReducedFIM(mt,list(freq=freq+df)))-
                  getMinEig(getReducedFIM(mt,list(freq=freq))))/df
  dfreq_exact = deig_dfreq(mt,freq) 
 
  expect_equal(as.numeric(dfreq_FD),as.numeric(dfreq_exact),tolerance =1e-4)
  
})