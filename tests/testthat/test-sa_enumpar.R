library(CVXR)
test_that("counting for n <= 5 works", {
  ref_val = c(1,1,2,3,5,7)
  expect_equal(sa_enumpar(5),ref_val)
})
