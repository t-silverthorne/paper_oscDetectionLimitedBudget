test_that("function evaluation and alpha dependence", {
  Nfine =2^6
  tau=seq(from=0,to=1,length.out=Nfine+1)
  tau=tau[1:(length(tau)-1)]
  xinds = sample(c(1:Nfine),10)
  tau[xinds]
  expect_gt(costfun_svdpower(tau[xinds],freqs=1,Amp=1,alpha=.05),
            costfun_svdpower(tau[xinds],freqs=1,Amp=1,alpha=.01))
})

test_that("takes min of freqs", {
  Nfine =2^6
  tau=seq(from=0,to=1,length.out=Nfine+1)
  tau=tau[1:(length(tau)-1)]
  xinds = sample(c(1:Nfine),10)
  tau[xinds]
  f1   = costfun_svdpower(tau[xinds],freqs=c(1),Amp=1,alpha=.05)
  f2   = costfun_svdpower(tau[xinds],freqs=c(2),Amp=1,alpha=.05)
  fall = costfun_svdpower(tau[xinds],freqs=c(1,2),Amp=1,alpha=.05)
  
  expect_equal(min(f1,f2),fall)
})