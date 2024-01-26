
test_that('exact power comparison',{
  t       = c(0,0.127,0.8,0.83)
  param   = list(Amp=3,freq=2.7,acro=pi)
  pwr     = eval_exact_power(t,param)
  pwr_ref =  0.155486657218106
  expect_equal(pwr,pwr_ref)
}) 


test_that('R power vs Monte Carlo',{
  Nmc=1e4
  t          = c(0:11)/11
  t          = t[1:10]
  param      = list(Amp=1.8,freq=2.3,acro=pi/3)
  pwr_exact  = eval_exact_power(t,param)
  pwr_mc=0
  nrep=1e1
  for (ii in c(1:nrep)){
    pwr_mc     = pwr_mc+eval_montecarlo_power(t,param,Nmc,alpha=.05)
  }
  pwr_mc=pwr_mc/nrep
  expect_equal(pwr_exact,pwr_mc,tolerance = 1e-2)
})