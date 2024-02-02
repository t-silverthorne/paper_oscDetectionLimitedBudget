test_that("function evaluation", {
  mt0=c(1:10)/10 - 1/10
  freqs=seq(from=1,to=24,length.out=24)
  Amp=1
  control=list(costfun_choice='svdpower',
               optim_method='L-BFGS-B',
               trace=0,#normally 6 if you want output
               REPORT=1,
               maxit=1)
  uu=opt_osc_power(mt0,freqs,Amp,control)
  expect_gt(uu$fvalue,0)
  expect_lt(uu$fvalue,1)
  
  control=list(costfun_choice='svdpower',
               optim_method='simul_anneal',
               tfun_choice = 'brownian-torus',
               tfun_mean   = 0, 
               tfun_sd     = .1, 
               trace=0,# make high if you want output
               REPORT=1,
               maxit=5)
  uu=opt_osc_power(mt0,freqs,Amp,control)
  expect_gt(uu$fvalue,0)
  expect_lt(uu$fvalue,1)
})
