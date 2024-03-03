require(devtools)
require(ggplot2)
require(annmatrix)
require(ggplotify)
require(patchwork)
require(parallel)
require(data.table)
require(stringr)
require(CVXR)
source('plotConfig.R')
load_all()

maxit=1e1
Nfreq=2^6
nrep = 10
mc_cores=10
control=list(
  N1min =  4,
  N1max = 20,
  maxit = maxit,
  trace = 1,
  REPORT = 50
)

nrep=1e1
freqs=seq(1,24,Nfreq)
mtdf=NULL
mtdf = c(1:nrep) %>% lapply(function(ind){
  control$N1min = floor(Nmeas/4)
  control$N1max = 3*floor(Nmeas/4)
  
  x0 =list(N1     = floor(Nmeas/4),
           N2     = 3*floor(Nmeas/4),
           shift2 = 1/4/Nmeas,
           scale2 = 1)
  res=solve_pin2lat(x0,freqs,control,cfuntype='ncp')
  data.frame(time=res$mtvalue,solver='pin2lat',Nmeas=Nmeas,
                                     ncp=-res$fvalue)
}) %>% rbindlist() %>% data.frame()
head(mtdf)
saveRDS(mtdf,'figs_for_paper/Fig1_aug_vs_cvxr_v2/solns_pin2lat.RDS')