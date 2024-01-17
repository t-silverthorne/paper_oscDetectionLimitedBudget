test_that("quadmats take correct value Linfty", {
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$costfun_type = 'Linfty'
  opts$Nfine = 5
  opts$Nfreq = 3
  opts$fmin  = 1
  opts$fmax  = 24
  Amat = make_quadmats(opts)
  
  # computed from matlab_checks/quadmats.m
  Aref1 = matrix(
   c(c(1.000000000000000  ,-0.809016994374949  , 0.309016994374952  ,  0.309016994374941  , -0.809016994374942),
   c(-0.809016994374949 ,  1.000000000000000 , -0.809016994374949 ,   0.309016994374952 ,   0.309016994374941),
   c(0.309016994374952  ,-0.809016994374949  , 1.000000000000000  , -0.809016994374949  ,  0.309016994374952),
   c(0.309016994374941  , 0.309016994374952  ,-0.809016994374949  ,  1.000000000000000  , -0.809016994374949),
   c(-0.809016994374942 ,  0.309016994374941 ,  0.309016994374952 ,  -0.809016994374949 ,   1.000000000000000)),
   nrow=opts$Nfine,byrow = T)
  Aref2 = replicate(25,{1}) %>% matrix(nrow=5)
  Aref3 = Aref1
  
  Alistref = list(Aref1,Aref2,Aref3)
  
  expect_equal(Amat,Alistref)
})

test_that("quadmats take correct value L1", {
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$costfun_type = 'L1'
  opts$Nfine = 8
  opts$Nfreq = 64 
  opts$fmin  = 1
  opts$fmax  = 24
  Amat = make_quadmats(opts)
  
  # computed from matlab_checks/quadmats.m
  Aref = matrix(
    c(c(1.000000000000000  ,-0.018683309633245  ,-0.000000000000000 ,  0.014535699883318  , 0.015625000000001   ,0.006734381679258  , 0.000000000000001 ,  0.004169471420661),
    c(-0.018683309633245 ,  1.000000000000000 , -0.018683309633245,  -0.000000000000000,   0.014535699883318  , 0.015625000000001,   0.006734381679259  , 0.000000000000001),
    c(-0.000000000000000 , -0.018683309633245 ,  1.000000000000000,  -0.018683309633244,  -0.000000000000000  , 0.014535699883318,   0.015625000000000  , 0.006734381679258),
    c(0.014535699883318  ,-0.000000000000000  ,-0.018683309633244 ,  1.000000000000000  ,-0.018683309633245  ,-0.000000000000000,   0.014535699883319 ,  0.015625000000001),
    c(0.015625000000001  , 0.014535699883318  ,-0.000000000000000 , -0.018683309633245  , 1.000000000000000  ,-0.018683309633245,  -0.000000000000001,   0.014535699883318),
    c(0.006734381679258  , 0.015625000000001  , 0.014535699883318 , -0.000000000000000  ,-0.018683309633245   ,1.000000000000000  ,-0.018683309633243  ,-0.000000000000000),
    c(0.000000000000001  , 0.006734381679259  , 0.015625000000000 ,  0.014535699883319  ,-0.000000000000001  ,-0.018683309633243,   1.000000000000000  ,-0.018683309633247),
    c(0.004169471420661  , 0.000000000000001  , 0.006734381679258 ,  0.015625000000001  , 0.014535699883318  ,-0.000000000000000,  -0.018683309633247,   1.000000000000000))
   ,nrow=opts$Nfine,byrow = T)
  
  expect_equal(Amat,Aref)
})

test_that('L1/Linfty cost function evaluation',{
  opts = make_default_opts(prob_size='small',
                           solver_type='simulanneal')
  opts$Nfine = 2^4
  opts$Nfreq = 2^6
  opts$Nmeas = 4
  opts$fmin  = 1
  opts$fmax  = 24
  x          = c(0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0)
  
  opts$costfun_type = 'L1'
  Amat       = make_quadmats(opts)
  L1_val     = sa_cfunpwr(x,Amat,opts) 
  L1_val2    = t(x)%*%Amat%*%x
  L1_ref     = 4.075165128996483 # from matlab_checks/cfun_eval.m
  expect_equal(L1_val,L1_val2)
  expect_equal(L1_val[1],L1_ref)
  
  opts$costfun_type = 'Linfty'
  Amat       = make_quadmats(opts)
  Linfty_val     = sa_cfunpwr(x,Amat,opts) 
  Linfty_ref     = 16 # from matlab_checks/cfun_eval.m
  expect_equal(Linfty_val,Linfty_ref)
  
})

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
  nrep=50
  for (ii in c(1:nrep)){
    pwr_mc     = pwr_mc+eval_montecarlo_power(t,param,Nmc,alpha=.05)
  }
  pwr_mc=pwr_mc/nrep
  expect_equal(pwr_exact,pwr_mc,tolerance = 1e-3)
})