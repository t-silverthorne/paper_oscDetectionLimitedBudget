test_that("compare with monte carlo", {
  param = list(freq=1+runif(1),
               Amp=1.2+.1*runif(1),
               acro=2*pi*runif(1))
  mt = c(1:15)/15-1/15
  
  al_val = sample(c(.01,.05),1) 
  pwr_exact = evalExactPower(mt,param,al_val)
  
  pwr_MC    = evalMonteCarloPower(mt,param,1e6,al_val)
  pwr_MC
  pwr_exact
  expect_equal(pwr_exact,pwr_MC,tolerance = 1e-2)
})

test_that('exact power comparison',{
  t       = c(0,0.127,0.8,0.83)
  param   = list(Amp=3,freq=2.7,acro=pi)
  pwr     = evalExactPower(t,param)
  pwr_ref =  0.155486657218106
  expect_equal(pwr,pwr_ref)
}) 

test_that("function evaluation, fast method matches slow", {
  param=list()
  param$Amp=2
  param$freq=2.1
  mt=c(0:9)/9
  mt=mt[1:(length(mt)-1)]
  Nacro=2^8
  acros = seq(from=0,to=2*pi,length.out=Nacro+1)
  acros = acros[1:Nacro]
  
  v1=min(apply(as.matrix(acros),1,function(acro){
    param$acro=acro
    return(evalExactPower(mt,param,alpha=.05))
    }))
  
  v2=evalWorstPower(mt,param$freq,param$Amp,alpha=.05)
  expect_equal(v1,v2,tolerance = 1e-3)
})