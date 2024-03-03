require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(CVXR)
load_all()
source('plotConfig.R')

tlim         = 5
Nmeasvals    = c(16,24)
threads_glob = 1

mtdf=NULL
control=list(
  Nmeas=24,
  Nfine=288,
  Nfreq=2^4,
  fmin=1,
  fmax=24,
  PreSolve=2,
  MIPFocus=3,
  cvxr_verbose=T,
  time_limit=tlim,
  maxit=1e9,
  MIPGapAbs=.01
)

for (Nmeas in Nmeasvals){
  control$Nmeas = Nmeas
  res=solve_cvxr(control,Threads=threads_glob,cfuntype='ncp')
  mtdf=rbind(mtdf,data.frame(time=res$mtvalue,solver='cvxr',Nmeas=Nmeas,
                             ncp=-res$fvalue))
}
head(mtdf)
saveRDS(mtdf,'figs_for_paper/Fig1_aug_vs_cvxr_v2/solns_cvxr.RDS')