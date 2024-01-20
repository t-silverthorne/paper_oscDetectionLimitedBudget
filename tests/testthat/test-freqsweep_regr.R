test_that("function evaluation", {
  mt=c(0:16)/16
  mt=mt[1:(length(mt)-1)]
  Nmc=1e1
  fmin=1
  fmax=2
  Amin=1
  Amax=2
  Xdat=make_simulated_data(mt,Nmc,Amin,Amax,fmin,fmax)
  Nfreq=3
  freqsweep_regr(mt,Xdat,fmin,fmax,Nfreq,return_type)
  #TODO: basic checks, amplitude in correct range, shape of data makes sense etc.
} 
)
