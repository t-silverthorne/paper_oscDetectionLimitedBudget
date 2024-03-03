require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(CVXR)
load_all()

Nmeas = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
maxit = 1e3
Nfreq = 2^10
nrep  =  10
mc_cores=Sys.getenv("SLURM_CPUS_PER_TASK")
control=list(
  N1min = NaN,
  N1max = NaN,
  maxit = maxit,
  trace = 1,
  REPORT = 50
)
control$N1min = floor(Nmeas/4)
control$N1max = 3*floor(Nmeas/4)
freqs=seq(1,24,length.out=Nfreq)
mtdf=NULL
mtdf = c(1:nrep) %>% lapply(function(ind){
  N1loc = sample(c(control$N1min,control$N1max),1) 
  x0 =list(N1     = N1loc,
           N2     = Nmeas - N1loc,
           shift2 = 1/2/N1loc,
           scale2 = 1)
  res=solve_pin2lat(x0,freqs,control,cfuntype='ncp')
  data.frame(time=res$mtvalue,solver='pin2lat',Nmeas=Nmeas,
                                     ncp=-res$fvalue)
}) %>% rbindlist() %>% data.frame()
head(mtdf)
saveRDS(mtdf,paste0('figs_for_paper/Fig1_aug_vs_cvxr_v2/solns_pin2lat_',Nmeas,'.RDS'))