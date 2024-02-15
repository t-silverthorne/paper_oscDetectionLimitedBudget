library(devtools)
library(CVXR)
library(gurobi)
load_all()

Nmeas         = 16 
Nfreq         = 2^9
Nfine         = 288
freqs_global  = seq(from=1,to=24,length.out=Nfreq)
trace_global  = 6
report_global = 1
maxit_global  = 3 
Amp_global    = 2

# control settings for each solver
ctrl_cts_arb_bfgs = list(costfun_choice='svdpower',
                      optim_method='L-BFGS-B',
                      trace=trace_global,#normally 6 if you want output
                      REPORT=report_global,
                      maxit=maxit_global)
var0_cts_arb_bfgs = c(1:Nmeas)/Nmeas - 1/Nmeas
resu_cts_arb_bfgs = opt_osc_power(dvar0   = var0_cts_arb_bfgs, 
                                  control = ctrl_cts_arb_bfgs,
                                  freqs   = freqs_global,
                                  Amp     = Amp_global)

ctrl_cts_arb_sa = list(costfun_choice='svdpower',
                      optim_method='simul_anneal',
                      tfun_choice='brownian-torus',
                      tfun_mean=0,
                      tfun_sd=.1,
                      trace=trace_global,#normally 6 if you want output
                      REPORT=report_global,
                      maxit=maxit_global)
var0_cts_arb_sa = c(1:Nmeas)/Nmeas - 1/Nmeas
resu_cts_arb_sa = opt_osc_power(dvar0   = var0_cts_arb_sa,
                                control = ctrl_cts_arb_sa,
                                freqs   = freqs_global,
                                Amp     = Amp_global)

#TODO: add wrapper to sweep over lattices
ctrl_cts_lat_bfgs = list(costfun_choice='svdpower_2lattice',
             optim_method='L-BFGS-B',
             trace=6,
             REPORT=1,
             maxit=1)
var0_cts_lat_bfgs = list(x0=c(shift2=0.5,scale1=.5,scale2=.5,shift1=0),
                         lat1 =NULL,
                         lat2 =NULL) 
resu_cts_lat_bfgs = wrapper_sweep_lattice(opt_osc_power,
                                Nvals    = c(4:6),
                                Nmeas    = Nmeas,
                                cts_flag = T,
                                dvar0    = var0_cts_lat_bfgs,
                                control  = ctrl_cts_lat_bfgs,
                                freqs    = freqs_global,
                                Amp      = Amp_global) 

ctrl_cts_lat_sa=list(costfun_choice = 'svdpower_2lattice',
             optim_method = 'simul_anneal',
             tfun_choice  = 'unif-with-bdry',
             tscale_unif_with_bdry=1/10,
             trace=5,
             REPORT=1,
             maxit=100)
var0_cts_lat_sa = list(x0=c(shift2=0.5,scale1=.5,scale2=.5,shift1=0),
                         lat1 =NULL,
                         lat2 =NULL) 
resu_cts_lat_sa = wrapper_sweep_lattice(opt_osc_power,
                                Nvals    = c(4:6),
                                Nmeas    = Nmeas,
                                cts_flag = T,
                                dvar0    = var0_cts_lat_sa,
                                control  = ctrl_cts_lat_sa,
                                freqs    = freqs_global,
                                Amp      = Amp_global) 

Nfine=288
ctrl_disc_arb_sa =list(costfun_choice='svdpower_discrete',
             optim_method='simul_anneal',
             tfun_choice='single-flip',
             trace=1,
             REPORT=1,
             maxit=50,
             fmin=1,fmax=24,Nfreq=2^10,
             Nfine=Nfine,Nmeas=Nmeas)
var0_disc_arb_sa = sample(c(1:Nfine),Nmeas) 
tau = c(1:Nfine)/Nfine -1/Nfine
resu_disc_arb_sa = opt_osc_power(dvar0  = var0_disc_arb_sa,
                                control = ctrl_disc_arb_sa,
                                freqs   = freqs_global,
                                Amp     = Amp_global,
                                tau     = tau)

ctrl_disc_arb_cvxr = list(costfun_choice='svdpower_discrete',
             optim_method='cvxr',
             trace=1,
             REPORT=1,
             maxit=1e9,time_limit=2,MIPGapAbs=.01,
             cvxr_verbose=T,
             costfun_type='Linfty',
             fmin=min(freqs_global),fmax=max(freqs_global),Nfreq=length(freqs_global),
             lattice_cstr='none',
             Nfine=Nfine,Nmeas=Nmeas)
resu_disc_arb_cvxr = opt_osc_power(dvar0  = NULL,
             control = ctrl_disc_arb_cvxr,
             freqs   = freqs_global,
             Amp     = Amp_global,
             tau     = tau)

#TODO: add wrapper to sweep over lattices
Nfine=288
ctrl_disc_lat_sa = list(costfun_choice = 'svdpower_2lattice_discrete',
             optim_method   = 'simul_anneal',
             tfun_choice    = 'unif-with-bdry-discrete',
             tscale         = 2,
             trace          = 6,
             REPORT         = 1,
             maxit          = 50)
var0_disc_lat_sa = list(x0 = c(dx1=1,dx2=1,xshift2=Nfine/2),
             N1 = NULL,
             N2 = NULL)
resu_disc_lat_sa = wrapper_sweep_lattice(opt_osc_power,
                                         Nvals    = c(4:6),
                                         Nmeas    = Nmeas,
                                         cts_flag = F,
                                         dvar0    = var0_disc_lat_sa,
                                         control  = ctrl_disc_lat_sa,
                                         freqs    = freqs_global,
                                         Amp      = Amp_global,
                                         tau      = tau)
