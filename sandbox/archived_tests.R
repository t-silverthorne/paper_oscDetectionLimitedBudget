#test_that("function evaluation", {
#  library(annmatrix)
#  library(dplyr)
#  library(data.table)
#  mt   = c(0:15)/15
#  mt   = mt[1:(length(mt)-1)]
#  Nmc  = 5e2
#  fmin = 1
#  fmax = 24 
#  Amin = 0.5 
#  Amax = 5 
#  Namp = 10 
#  true_freq_vals  = c(1,2,12) %>% as.list()
#  Ampvals         = seq(from=Amin,to=Amax,length.out=Namp) %>% as.list()
#  Nfreq_regr_vals = c(10,15,20) %>% as.list()
#  start=Sys.time()
#  uu=benchmark_fdr(mt,Nmc,fmin,fmax,Ampvals,Nfreq_regr_vals,'fdr')
#  end=Sys.time()
#  end-start
#  library(ggplot2)
#  library(patchwork)
#  uu %>% head()
#  
#  #p1=uu %>% ggplot(aes(x=Amp,y=pdetect_p,color=as.factor(true_freq),group=true_freq))+
#  #  geom_line()+facet_wrap(~Nfreq_regr,nrow=3) +theme(legend.position='none')
#  #p2=uu %>% ggplot(aes(x=Amp,y=pdetect_q,color=as.factor(true_freq),group=true_freq))+
#  #  geom_line()+facet_wrap(~Nfreq_regr,nrow=3) +theme(legend.position='none')
#  #
#  #p1+p2
#})

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
  nrep=1e1
  for (ii in c(1:nrep)){
    pwr_mc     = pwr_mc+eval_montecarlo_power(t,param,Nmc,alpha=.05)
  }
  pwr_mc=pwr_mc/nrep
  expect_equal(pwr_exact,pwr_mc,tolerance = 1e-2)
})test_that("compare to ref", {
lat1=c(1:4)/4-1/4
lat2=c(1:4)/4-1/4

lat_ref = c(1:8)/8-1/8
expect_equal(sort(convert_2lattice_to_state(0,1/8,1,1,lat1,lat2)),lat_ref)
})

test_that("scaling shifting goes forward", {
lat1=c(1:4)/4-1/4
lat2=c(1:4)/4-1/4

lat_ref = c(1:8)/8-1/8
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(1/128,0,1,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,1/128,1,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,0,1.01,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,0,1,1.01,lat1,lat2)),)
})test_that("function evaluation, alpha handling", {
  N1=4
  N2=4
  dx1=1
  dx2=4
  xshift2=60
  Nfine=1e2
  xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
  e1=costfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,c(1:Nfine),c(1,5),
                                     Amp=1,alpha=.05)
  e2=costfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,c(1:Nfine),c(1,1,1),
                                     Amp=1,alpha=.01)
  expect_gt(e1,e2)
})
test_that("function evaluation, alpha handling", {
shift1=0
shift2=.1
scale1=1
scale2=1
lat1=c(1:5)/5-1/5
lat2=c(1:5)/5-1/5
expect_gte(
  costfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,c(1,2,3),.05),
  costfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,c(1,2,3),.01))
})

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
    return(eval_exact_power(mt,param,alpha=.05))
    }))
  
  v2=costfun_svdpower(mt,param$freq,param$Amp,alpha=.05)
  expect_equal(v1,v2,tolerance = 1e-3)
})
test_that("function evaluation", {
  mt = runif(30)
  expect_no_error(deig_dfreq(mt,runif(1)))
})

test_that('agrees with finite difference',{
  set.seed(1)
  mt = runif(20)
  
  freq = 2.5
  df   = 1e-5
  
  dfreq_FD    = (getMinEig(getReducedFIM(mt,list(freq=freq+df)))-
                  getMinEig(getReducedFIM(mt,list(freq=freq))))/df
  dfreq_exact = deig_dfreq(mt,freq) 
 
  expect_equal(as.numeric(dfreq_FD),as.numeric(dfreq_exact),tolerance =1e-4)
  
})test_that("multiplication works", {
  mt=runif(16)
  
  # cost function value without any regularization
  c1=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),cfuntype = 'ncp')
  c1b=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=0,cfuntype = 'ncp')
  c1c=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regFder=0,cfuntype = 'ncp')
  c1d=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=0,regFder=0,cfuntype = 'ncp')
  expect_equal(c1,c1b)
  expect_equal(c1b,c1c)
  expect_equal(c1c,c1d)
 
  # add mean freq regularization 
  c2=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=1,cfuntype = 'ncp')
  expect_gt(c2,c1)
 
  # add frequency derivative regularization 
  c3=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=1,regFder=1,cfuntype = 'ncp')
  expect_lt(c3,c2)
  
  expect_error(costfun_svdpower(mt,freqs=c(1,2,3),regFder=10),'Use cfuntype')
})
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
test_that("function evaluation", {
  Nfine=2^7
  tau=c(1:Nfine)/Nfine - 1/Nfine 
  freqs=seq(from=1,to=24,length.out=24)

  dvar0=list(x0 = c(dx1=1,dx2=1,xshift2=Nfine/2),
  N1      = 8,
  N2      = 6)

  control=list(costfun_choice = 'svdpower_2lattice_discrete',
               optim_method   = 'simul_anneal',
               tfun_choice    = 'unif-with-bdry-discrete',
               tscale         = 2,
               trace          = 0,
               REPORT         = 1,
               maxit          = 50)
  xout=opt_osc_power(dvar0=dvar0,freqs=freqs,control=control,tau=tau)

  expect_gt(-xout$fvalue,0)
  expect_lt(-xout$fvalue,1)
})
test_that("costfun evaluates", {
  N1=40
  N2=4
  scale2=.5
  shift2=1/20
  freqs=seq(from=1,to=24,length.out=2^8)
  r1=costfun_auglattice(N1,N2,shift2,scale2,freqs,cfuntype='power',Amp=1.5)
  r2=costfun_auglattice(N1,N2,shift2,scale2,freqs,cfuntype='ncp')
  expect_lt(r1,r2)
})

test_that("transition function conserves Nmeas", {
  Nmeas=24
  N1=12
  N2=12
  scale2=.5
  shift2=1/20
  control=list(N1min=floor(Nmeas/4),N1max=Nmeas)
  for (ii in c(1:100)){
    xnew=tfun_auglattice(N1,N2,shift2,scale2,control=control)
    N1=xnew[1]
    N2=xnew[2]
    shift2=xnew[3]
    scale2=xnew[4]
  }
  expect_equal(N1+N2,Nmeas)
})

test_that("solution function evaluates", {
  N1=12
  dvar0=list(N1=12,
             N2=12,
             shift2=1/2/N1,
             scale2=1)
  freqs=seq(from=1,to=24,length.out=2^8)
  control=list(N1min=8,N1max=16,trace=0,REPORT=0,maxit=1e2)
  xout=solve_auglattice(dvar0=dvar0,
                        freqs=freqs,
                        control=control,
                        cfuntype='ncp')
  expect_lt(xout$value,0)
  expect_equal(xout$par[1]+xout$par[2],24)
})

test_that('opt_osc_power wrapper works on auglattice',{
  N1=16
  dvar0=list(N1=N1,
             N2=8,
             shift2=1/2/N1,
             scale2=1)
  freqs=seq(from=1,to=24,length.out=2^8)
  control=list(N1min=8,N1max=16,trace=0,REPORT=0,maxit=1e1,
               costfun_choice='auglattice',cfuntype='ncp')
  expect_no_error(opt_osc_power(dvar0=dvar0,freqs=freqs,control=control) )
})

test_that('second test of evaluation',{
  gset = list(
    maxit_sa          = 1e3,     # max iterations for simulated annealing 
    timelimit_cvxr    = 5*60,    # 20 seconds of compute time
    Nmeas             = NaN,     # measurement budget
    Nfreq             = 2^3,     # number of frequencies to use
    Nfine             = 288,     # corresponds to 5 minute sampling
    trace_global      = 6,       # controls how often to report 
    report_global     = 10,      # how much reporting you want
    Amp_global        = NaN,     # amp is irrelevant, not used here
    Nmin_2lat         = NaN,     # will be updated in inner loop
    Nmax_2lat         = NaN      # will be updated in inner loop
  )
  freqs=seq(from=1,to=24,length.out=gset$Nfreq)
  Nmeas=16
  gset$Nmeas=Nmeas
  gset$Nmin_2lat=floor(Nmeas/4)
  gset$Nmax_2lat=3*floor(Nmeas/4)
  dvar0=list(N1=floor(Nmeas/4),
             N2=3*floor(Nmeas/4),
             shift2=1/4/Nmeas,
             scale2=1)
  control=list(N1min=gset$Nmin_2lat,N1max=gset$Nmax_2lat,
               trace=1,REPORT=1,maxit=gset$maxit_sa,
               costfun_choice='auglattice')
  #res=opt_osc_power(dvar0=dvar0,freqs=freqs,control=control,
  #           cfuntype='ncp',
  #           regFder=10)
  expect_no_error(opt_osc_power(dvar0=dvar0,freqs=freqs,control=control,
                    cfuntype='ncp'))
})test_that("function evaluation", {
  shift1=0
  shift2=runif(1)*.5
  scale1=.5
  scale2=.5
  lat1=c(1:5)/5 - 1/5
  lat2=c(1:10)/10 -1/10
  Amp=1
  freqs=seq(from=1,to=24,length.out=24)
  control=list(costfun_choice='svdpower_2lattice',
               optim_method='L-BFGS-B',
               trace=6,
               REPORT=1,
               maxit=1)
  x0=c(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2)
  xout=opt_osc_power(dvar0=list(x0=x0,lat1=lat1,lat2=lat2),freqs=freqs,control=control,Amp=Amp)
  expect_gt(-xout$fvalue,0)
  expect_lt(-xout$fvalue,1)
 
  control$optim_method ='simul_anneal' 
  control$maxit        = 50
  control$tfun_choice  = 'unif-with-bdry'
  control$tscale_unif_with_bdry=1/10
  x0=c(shift1=shift1,shift2=shift2,scale1=scale1,scale2=scale2)
  xout=opt_osc_power(dvar0=list(x0=x0,lat1=lat1,lat2=lat2),freqs=freqs,control=control,Amp=Amp)
  expect_gt(-xout$fvalue,0)
  expect_lt(-xout$fvalue,1)
  
})
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
test_that("function evaluation", {
  mt0=c(1:10)/10 - 1/10
  freqs=seq(from=1,to=24,length.out=24)
  Amp=1
  control=list(costfun_choice='svdpower',
               optim_method='L-BFGS-B',
               trace=0,#normally 6 if you want output
               REPORT=1,
               maxit=1)
  uu=opt_osc_power(dvar0=mt0,freqs=freqs,control=control)
  expect_gt(-uu$fvalue,0)
  expect_lt(-uu$fvalue,1)
  
  control=list(costfun_choice='svdpower',
               optim_method='simul_anneal',
               tfun_choice = 'brownian-torus',
               tfun_mean   = 0, 
               tfun_sd     = .1, 
               trace=0,# make high if you want output
               REPORT=1,
               maxit=5)
  uu=opt_osc_power(dvar0=mt0,freqs=freqs,control=control)
  expect_gt(-uu$fvalue,0)
  expect_lt(-uu$fvalue,1)
})
test_that("generates valid state", {
N1=sample(c(8:12),1) 
N2=sample(c(6:12),1)
dx1=1
dx2=1
xshift2=50
Nfine=200
control=list(tfun_choice='unif-with-bdry-discrete')
control$tscale=2
expect_false(
  as.logical(anyDuplicated(xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)))
  )
for (ii in c(1:1e3)){
  L=tfun_2lattice_svdpower_discrete(N1,dx1,N2,dx2,xshift2,Nfine,control)
  dx1=L$dx1
  dx2=L$dx2
  xshift2=L$xshift2
}
xinds=xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
expect_false(
  as.logical(anyDuplicated(xinds))
  )
expect_true(all(xinds>0))
expect_true(all(xinds<=Nfine))
})
test_that("points stay inside lattice", {
  control=list(tfun_choice='unif-with-bdry')
  control$tscale_unif_with_bdry=1/10
  shift1=runif(1)/12
  shift2=runif(1)/12
  scale1=runif(1)
  scale2=runif(1)
  N1=sample(c(8:12),1)
  N2=sample(c(8:12),1)
  lat1= c(1:N1)/N1 -1/N1
  lat2= c(1:N2)/N2 -1/N2
  for (ii in c(1:100)){
    L=tfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,control)
    scale1=L$scale1
    scale2=L$scale2
    shift2=L$shift2
  }
  mt=convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2)
  expect_false(any(mt>1))
  expect_false(any(mt<0))
})
test_that("function conserves inds, moves at least one", {
  Nfine = sample(c(50:100),1)
  Nmeas = sample(c(1,20),1)
  tau      = c(1:Nfine)/Nfine -1/Nfine
  xinds    = sample(c(1:Nfine),Nmeas)
  control  = list(tfun_choice='single-flip')
  expect_equal(length(unique(tfun_svdpower_discrete(xinds,Nfine,control))),Nmeas)
  expect_equal(length(which(abs(xinds-tfun_svdpower_discrete(xinds,Nfine,control))>0)),1)
})
test_that("basic reference test", {
N1=5
N2=10
dx1=1
dx2=2
xshift2=6

ref=sort(c(1+dx1*c(0:(N1-1)),1+xshift2+dx2*c(0:(N2-1))))
xinds=xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2)
expect_equal(xinds,ref)
})
