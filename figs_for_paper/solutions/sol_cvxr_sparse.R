require(devtools)
require(annmatrix)
require(stringr)
require(CVXR)
require(gurobi)
load_all()

Nfreq        =  47
tlim         =  60*60 
Nmeas=as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
threads_glob = Sys.getenv("SLURM_CPUS_PER_TASK") 

mtdf=NULL
control=list(
  Nmeas=Nmeas,
  Nfine=288,
  Nfreq=Nfreq,
  fmin=1,
  fmax=24,
  PreSolve=-1,
  MIPFocus=3,
  cvxr_verbose=T,
  time_limit=tlim,
  maxit=1e9,
  MIPGapAbs=.01,
  NodefileStart=0.45,
  prob_formulation='max_elemwise'
  )

res=solve_cvxr(control,Threads=threads_glob,cfuntype='ncp')
mtdf=rbind(mtdf,data.frame(time=res$mtvalue,solver='cvxr',Nmeas=Nmeas,
                            ncp=-res$fvalue))
head(mtdf)
saveRDS(mtdf,paste0('figs_for_paper/solutions/solns_cvxr_sparse_',Nmeas,'.RDS'))
