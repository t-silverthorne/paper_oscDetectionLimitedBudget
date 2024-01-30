test_that("function eval", {
  Nfine=2^7
  tau=c(1:Nfine)/Nfine - 1/Nfine 
  xinds = sample(c(1:Nfine),20) 
  
  freqs=seq(from=1,to=24,length.out=24)
  Amp=1
  control=list(costfun_choice='svdpower_discrete',
               optim_method='simul_anneal',
               tfun_choice ='single-flip',
               trace=1,
               REPORT=1,
               maxit=1)
  xout=opt_osc_power(mt0=xinds,freqs=freqs,Amp=Amp,control=control,tau=tau)
  expect_gt(xout$value,0)
  expect_lt(xout$value,1)
})
