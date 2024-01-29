test_that("function evaluation, alpha handling", {
shift1=0
shift2=.1
scale1=1
scale2=1
lat1=c(1:5)/5-1/5
lat2=c(1:5)/5-1/5
expect_gte(
  costfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,c(1,2,3),.05),
  costfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,c(1,2,3),.01))
})

