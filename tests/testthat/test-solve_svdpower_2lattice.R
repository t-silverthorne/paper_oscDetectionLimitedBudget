test_that("function evaluation", {
  shift1=0
  shift2=runif(1)*.5
  scale1=.5
  scale2=.5
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
  expect_gt(xout$fvalue,0)
  expect_lt(xout$fvalue,1)
 
  control$optim_method ='simul_anneal' 
  control$maxit        = 50
  control$tfun_choice  = 'unif-with-bdry'
  control$tscale_unif_with_bdry=1/10
  x0=c(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2)
  xout=opt_osc_power(dvar0=list(x0=x0,lat1=lat1,lat2=lat2),freqs,Amp,control)
  expect_gt(xout$fvalue,0)
  expect_lt(xout$fvalue,1)
  
})
