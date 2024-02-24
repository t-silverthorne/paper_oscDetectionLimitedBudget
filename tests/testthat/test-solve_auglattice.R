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
})