test_that("function evaluates", {
  fmin  = 1
  fmax  = 3.8 
  Nf    = 2
  Nacro = 2^6
  Nmeas = 12 
  Amp   = 1.6
  ctrl = list(verbose=T,smooth=F)
  uu=wrap_optimsa(Nmeas,fmin,fmax,Nf,Nacro,Amp)
  par(mfrow=c(1,2))
  plot(-uu$Sout,type='l')
  plot(uu$Tout,type='l')
  expect_gte(-uu$Sout[length(uu$Sout)],0)
})
