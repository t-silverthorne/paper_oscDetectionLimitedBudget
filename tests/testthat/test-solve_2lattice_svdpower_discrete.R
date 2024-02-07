test_that("function evaluation", {
  Nfine=2^7
  tau=c(1:Nfine)/Nfine - 1/Nfine 
  freqs=seq(from=1,to=24,length.out=24)

  dvar0=list(x0 = c(dx1=1,dx2=1,xshift2=Nfine/2),
  N1      = 8,
  N2      = 6)

  control=list(costfun_choice = 'svdpower_2lattice_discrete',
               optim_method   = 'simul_anneal',
               tfun_choice    = 'unif-with-bdry-discrete',
               tscale         = 2,
               trace          = 0,
               REPORT         = 1,
               maxit          = 50)
  xout=opt_osc_power(dvar0=dvar0,freqs=freqs,Amp=2,control=control,tau=tau)

  expect_gt(xout$fvalue,0)
  expect_lt(xout$fvalue,1)
})
