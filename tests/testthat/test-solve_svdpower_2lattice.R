test_that("function evaluation", {
  shift1=0
  shift2=runif(1) 
  scale1=1
  scale2=1
  lat1=c(1:5)/5 - 1/5
  lat2=c(1:10)/10 -1/10
  Amp=1
  freqs=seq(from=1,to=24,length.out=24)
  control=list(costfun_choice='svdpower_2lattice',
               optim_method='L-BFGS-B',
               trace=6,
               REPORT=1,
               maxit=1)
  x0=c(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2)
  xout=opt_osc_power(dvar0=list(x0=x0,lat1=lat1,lat2=lat2),freqs,Amp,control)
  expect_gt(xout$value,0)
  expect_lt(xout$value,1)
 
  control$optim_method ='simul_anneal' 
  control$maxit        = 5
  control$tfun_choice  = 'default'
  control$mean_scale1  = 0 
  control$mean_scale2  = 0 
  control$mean_shift2  = 0 
  control$sd_scale1    = .1
  control$sd_scale2    = .1
  control$sd_shift2    = .1
  x0=c(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2)
  xout=opt_osc_power(dvar0=list(x0=x0,lat1=lat1,lat2=lat2),freqs,Amp,control)
  expect_gt(xout$value,0)
  expect_lt(xout$value,1)
  
})
