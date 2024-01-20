#test_that("function evaluation", {
#  load_all('.')
#  mt=c(0:16)/16
#  mt=mt[1:(length(mt)-1)]
#  Nmc=1e1
#  fmin=1
#  fmax=2
#  Amin=1
#  Amax=2
#  Xdat=make_simulated_data(mt,Nmc,Amin,Amax,fmin,fmax)
#  Nfreq=3
#  L        = freqsweep_regr(mt,Xdat,fmin,fmax,Nfreq,return_type)
#  qvals    = matrix_1d_padjust(L$pvalues,1,'fdr')
#  dom_freq = extract_dominant_freq(qvals,L$amps,L$acros)
#  #TODO: basic checks, amplitude in correct range, shape of data makes sense etc.
#} 
#)
