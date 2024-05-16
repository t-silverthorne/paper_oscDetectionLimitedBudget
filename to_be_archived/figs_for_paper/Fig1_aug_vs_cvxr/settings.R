require(devtools)
require(ggplot2)
require(annmatrix)
require(ggplotify)
require(patchwork)
require(parallel)
require(data.table)
require(stringr)
require(CVXR)
load_all()
gset = list(
  maxit_sa          = 1e3,     # max iterations for simulated annealing 
  timelimit_cvxr    = 5*60,    # units of seconds of compute time
  Nmeas             = NaN,     # measurement budget
  Nfreq             = 2^10,     # number of frequencies to use
  Nfine             = 288,     # corresponds to 5 minute sampling
  trace_global      = 6,       # controls how often to report 
  report_global     = 10,      # how much reporting you want
  Amp_global        = NaN,     # amp is irrelevant, not used here
  Nmin_2lat         = NaN,     # will be updated in inner loop
  Nmax_2lat         = NaN      # will be updated in inner loop
)

fmin=1
fmax=24
freqs=seq(from=1,to=24,length.out=gset$Nfreq)
Nmeasvals = c(16,24,32)

N_amp_plt = 50  # 50
N_per_plt = 100 # 1e2

NmcFDR    = 5e2 # 5e2
NacroFDR  = 2^5 # 2^5
NfreqFDR  = 2^5 # 2^5