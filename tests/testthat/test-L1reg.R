require(devtools)
require(CVXR)
test_that("L1 regularization has effect on continuous", {
  mt0=c(1:10)/10 - 1/10
  freqs=seq(from=1,to=24,length.out=24)
  Amp=1
  control=list(costfun_choice='svdpower',
               optim_method='L-BFGS-B',
               trace=0,#normally 6 if you want output
               REPORT=1,
               maxit=5)
  set.seed(1)
  uu=opt_osc_power(mt0,freqs,control,regL1=0,Amp=Amp)
  set.seed(1)
  vv=opt_osc_power(mt0,freqs,control,regL1=0,Amp=Amp)
  set.seed(1)
  ww=opt_osc_power(mt0,freqs,control,regL1=10,Amp=Amp)
  
  expect_equal(uu$mtvalue,vv$mtvalue)
  expect_equal(sum((uu$mtvalue-vv$mtvalue)^2),0)
  expect_gt(sum((uu$mtvalue-ww$mtvalue)^2),0)
  
})

#TODO: fix CVXR L1 regularization
#test_that('L1 regularization has effect on CVXR',{
#  Nmeas=20
#  Nfine = 2^5 
#  freqs=seq(from=1,to=24,length.out=24)
#  Amp_global = 2
#  freqs_global  = seq(from=1,to=24,length.out=8)
#  tau = c(1:Nfine)/Nfine -1/Nfine
#  ctrl_disc_arb_cvxr = list(costfun_choice='svdpower_discrete',
#              optim_method='cvxr',
#              maxit=1e9,time_limit=2,MIPGapAbs=.01,
#              cvxr_verbose=F,costfun_type='Linfty',
#              fmin=min(freqs),
#              fmax=max(freqs),Nfreq=length(freqs),
#              lattice_cstr='none',Nfine=Nfine,Nmeas=Nmeas)
#  res1=opt_osc_power(control = ctrl_disc_arb_cvxr,
#               freqs   = freqs_global,
#               Amp     = Amp_global,
#               tau     = tau)
#  res2=opt_osc_power(control = ctrl_disc_arb_cvxr,
#               freqs   = freqs_global,
#               Amp     = Amp_global,
#               tau     = tau,
#               regL1   = 1e1)
#  expect_gt(sum((res1$mtvalue-res2$mtvalue)^2),0)
#})
