test_that("function eval", {
  Nfine=2^7
  tau=c(1:Nfine)/Nfine - 1/Nfine 
  xinds = sample(c(1:Nfine),20) 
  
  freqs=seq(from=1,to=24,length.out=24)
  Amp=1
  control=list(costfun_choice='svdpower_discrete',
               optim_method='simul_anneal',
               tfun_choice ='single-flip',
               Nfine=length(tau),
               trace=1,
               REPORT=1,
               maxit=2)
  xout=opt_osc_power(dvar0=xinds,freqs=freqs,control=control,tau=tau,Amp=Amp)
  expect_gt(-xout$fvalue,0)
  expect_lt(-xout$fvalue,1)
  
  library(CVXR)
  library(gurobi)
  control=list(costfun_choice='svdpower_discrete',
               optim_method='cvxr',
               trace=1,
               REPORT=1,
               maxit=1e9,time_limit=2,MIPGapAbs=.01,
               cvxr_verbose=F,
               costfun_type='Linfty',
               fmin=1,fmax=24,Nfreq=8,
               lattice_cstr='none',
               Nfine=2^4,Nmeas=16)#TODO: find out why Nfine=2^5 gives NA
  xout=opt_osc_power(dvar0=xinds,freqs=freqs,Amp=Amp,control=control,tau=tau)
  expect_gt(-xout$fvalue,0)
  expect_lt(-xout$fvalue,1)
  
})
