require(devtools)
require(annmatrix)
require(stringr)
require(CVXR)
require(gurobi)
load_all()

Nfreq        =  47
tlim         =  2*60*60 
Nmeas=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
threads_glob = Sys.getenv("SLURM_CPUS_PER_TASK") 

mtdf=NULL
control=list(
  Nmeas=Nmeas,
  Nfine=144,
  Nfreq=Nfreq,
  fmin=1,
  fmax=24,
  PreSolve=2,
  MIPFocus=3,
  cvxr_verbose=T,
  time_limit=tlim,
  maxit=1e9,
  MIPGapAbs=1e-12,
  MIPGap=1e-5,
  NodefileStart=Inf
  )

#res=solve_cvxr_spt(control,Threads=threads_glob,
#                   use_spt=F,
#                   drts=NULL,
#                   pads=NULL)
#saveRDS(res,paste0('figs_for_paper/support_constrained_solutions/solns_cvxr_',Nmeas,'.RDS'))

drts=6*2
pads=0
res=solve_cvxr_spt(control,Threads=threads_glob,
                   use_spt=T,
                   drts=drts,
                   pads=pads)
saveRDS(res,paste0('figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_',drts,
                   '_pads_',pads,
                   '_Nmeas_',Nmeas,'.RDS'))



