test_that("takes min of freqs", {
  Nfine =2^6
  tau=seq(from=0,to=1,length.out=Nfine+1)
  tau=tau[1:(length(tau)-1)]
  xinds = sample(c(1:Nfine),10)
  tau[xinds]
  f1   = evalWorstPower(tau[xinds],freqs=c(1),Amp=2)
  f2   = evalWorstPower(tau[xinds],freqs=c(2),Amp=2)
  fall = evalWorstPower(tau[xinds],freqs=c(1,2),Amp=2)
  expect_equal(min(f1,f2),fall)
})
