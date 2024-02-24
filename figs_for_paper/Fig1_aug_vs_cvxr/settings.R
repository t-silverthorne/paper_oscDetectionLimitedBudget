gset = list(
  maxit_sa          = 1e3,     # max iterations for simulated annealing 
  timelimit_cvxr    = 5*60,    # 20 seconds of compute time
  Nmeas             = NaN,     # measurement budget
  Nfreq             = 2^8,     # number of frequencies to use
  Nfine             = 288,     # corresponds to 5 minute sampling
  trace_global      = 6,       # controls how often to report 
  report_global     = 10,      # how much reporting you want
  Amp_global        = NaN,     # amp is irrelevant, not used here
  Nmin_2lat         = NaN,     # will be updated in inner loop
  Nmax_2lat         = NaN      # will be updated in inner loop
)
