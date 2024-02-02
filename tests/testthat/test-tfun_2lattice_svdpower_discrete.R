test_that("generates valid state", {
N1=sample(c(8:12),1) 
N2=sample(c(6:12),1)
dx1=1
dx2=1
xshift2=50
Nfine=200
control=list(tfun_choice='unif-with-bdry-discrete')
control$tscale=2
expect_false(
  as.logical(anyDuplicated(xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)))
  )
for (ii in c(1:1e3)){
  L=tfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,Nfine,control)
  dx1=L$dx1
  dx2=L$dx2
  xshift2=L$xshift2
}
xinds=xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
expect_false(
  as.logical(anyDuplicated(xinds))
  )
expect_true(all(xinds>0))
expect_true(all(xinds<=Nfine))
})
