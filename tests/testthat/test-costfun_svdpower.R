test_that("fast method matches slow", {
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
    return(eval_exact_power(mt,param,alpha))
    }))
  
  v2=costfun_svdpower(mt,param$freq,param$Amp,alpha=.05)
  expect_equal(v1,v2,tolerance = 1e-3)
})
