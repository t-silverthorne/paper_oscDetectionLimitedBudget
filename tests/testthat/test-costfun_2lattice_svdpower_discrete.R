test_that("function evaluation, alpha handling", {
  N1=4
  N2=4
  dx1=3
  dx2=4
  xshift2=1
  Nfine=1e2
  xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
  e1=costfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,c(1:Nfine),c(1,5),
                                     Amp=1,alpha=.05)
  e2=costfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,c(1:Nfine),c(1,1,1),
                                     Amp=1,alpha=.01)
  expect_gt(e1,e2)
})
