require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(CVXR)
require(gurobi)
load_all()

Nfreq        = 2^8
tlim         =  10
Nmeas=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
threads_glob = Sys.getenv("SLURM_CPUS_PER_TASK") 

mtdf=NULL
control=list(
  Nmeas=Nmeas,
  Nfine=288,
  Nfreq=Nfreq,
  fmin=1,
  fmax=24,
  PreSolve=2,
  MIPFocus=3,
  cvxr_verbose=T,
  time_limit=tlim,
  maxit=1e9,
  MIPGapAbs=.01
)

res=solve_cvxr(control,Threads=threads_glob,cfuntype='ncp')
mtdf=rbind(mtdf,data.frame(time=res$mtvalue,solver='cvxr',Nmeas=Nmeas,
                            ncp=-res$fvalue))
head(mtdf)
saveRDS(mtdf,paste0('figs_for_paper/Fig1_aug_vs_cvxr_v2/solns_cvxr_',Nmeas,'.RDS'))