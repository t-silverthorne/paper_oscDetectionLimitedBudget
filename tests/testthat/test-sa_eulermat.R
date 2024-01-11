library(CVXR)
test_that("naive test: Emat[:,1]*Emat[:,2] <= m", {
  m = sample(c(5:10),1) 
  Emat = sa_eulermat(m)
  counts = Emat %>% lapply(function(x){x[1]*x[2]}) %>% unlist()
  expect_lte(max(counts),m)
  expect_gte(min(counts),1)
})

#TODO: check probabilities Emat[:,3] are correct by doing small example by hand