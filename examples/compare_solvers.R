library(devtools)
library(CVXR)
library(gurobi)
load_all()

Nmeas       = 16 
Nfreq       = 2^10
Nfine       = 288
freqs       = seq(from=0,to=24,length.out=Nfreq)
trace_glob  = 6
report_glob = 1
maxit_glob  = 3 
Amp_global  = 2

# control settings for each solver
ctrl_cts_arb_bfgs = list(costfun_choice='svdpower',
                      optim_method='L-BFGS-B',
                      trace=trace_glob,#normally 6 if you want output
                      REPORT=report_glob,
                      maxit=maxit_glob)
var0_cts_arb_bfgs = c(1:Nmeas)/Nmeas - 1/Nmeas
resu_cts_arb_bfgs = opt_osc_power(dvar0   = var0_cts_arb_bfgs, 
                                  control = ctrl_cts_arb_bfgs,
                                  freqs   = freqs,
                                  Amp     = Amp_global)

ctrl_cts_arb_sa = list(costfun_choice='svdpower',
                      optim_method='simul_anneal',
                      tfun_choice='brownian-torus',
                      tfun_mean=0,
                      tfun_sd=.1,
                      trace=trace_glob,#normally 6 if you want output
                      REPORT=report_glob,
                      maxit=maxit_glob)
var0_cts_arb_sa = c(1:Nmeas)/Nmeas - 1/Nmeas
resu_cts_arb_sa = opt_osc_power(dvar0   = var0_cts_arb_sa,
                                control = ctrl_cts_arb_sa,
                                freqs   = freqs,
                                Amp     = Amp_global)

#TODO: add wrapper to sweep over lattices
ctrl_cts_lat_bfgs 
var0_cts_lat_bfgs

ctrl_cts_lat_sa
var0_cts_lat_sa

ctrl_disc_arb_sa
var0_disc_arb_sa

ctrl_disc_arb_cvxr
var0_disc_arb_cvxr

#TODO: add wrapper to sweep over lattices
ctrl_disc_lat_sa
var0_disc_lat_sa

# run all solvers



# compare results