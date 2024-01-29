test_that("basic reference test", {
N1=5
N2=10
dx1=1
dx2=2
xshift2=6

ref=sort(c(1+dx1*c(0:(N1-1)),1+xshift2+dx2*c(0:(N2-1))))
xinds=xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
expect_equal(xinds,ref)
})
