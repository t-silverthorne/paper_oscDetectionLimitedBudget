require(matrixTests)
test_that("compare with monte carlo", {
  param = list(freq=1+runif(1),
               Amp=1+.1*runif(1),
               acro=2*pi*runif(1))
  mt = c(1:25)/25-1/25
  Nmc   = 5e3
 
  # optimal choice of Nperm based on Boos Zhang heuristic 
  malpha  = 10
  al_val  = 1/malpha
  Npcands = malpha*c(1:1000)-1
  Nperm   = Npcands[which.min(abs(Npcands-8*sqrt(Nmc)))]
  
  pwr_exact = evalExactPower(mt,param,al_val)
  pwr_MC1 = evalMonteCarloPower(mt,param,Nmc,al_val,method='Ftest')
  pwr_MC2 = evalMonteCarloPower(mt,param,Nmc,al_val,method='perm',Nperm=Nperm)
  expect_equal(pwr_exact,pwr_MC1,tolerance = 1e-2)
  expect_equal(pwr_exact,pwr_MC2,tolerance = 1e-2)
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