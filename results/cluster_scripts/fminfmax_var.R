require(dplyr)
require(gurobi)
source('results/source/powerChord.R')

# configure settings grid
out_loc_glob='results/data/fminfmax/'
Nmeas_vals = c(16,24,32,48)
fmin_vals  = c(1,2,4,6,8,10,12) 
df_vals    = seq(1,23,3) 
setgs      = expand.grid(Nmeas = Nmeas_vals,
                         fmin  = fmin_vals,
                         df    = df_vals)

# only consider priors with at most 1hr rhythms
setgs      = setgs[(setgs$fmin+setgs$df)<=24,] 

ii           = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
num_threads  = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

sol=powerChord(Nmeas       = setgs[ii,]$Nmeas,
               fmin        = setgs[ii,]$fmin,
               fmax        = setgs[ii,]$fmin+setgs[ii,]$df,
               Nfreq       = 256,
               drts        = Inf,
               w_reg       = 0,
               num_threads = num_threads,
               tlim        = 60*60)

saveRDS(sol,file=paste0(out_loc_glob,'sol_gur_',
                        '_Nmeas_',setgs[ii,]$Nmeas,
                        '_fmin_',setgs[ii,]$fmin,
                        '_df_',setgs[ii,]$df,                     
                        '.RDS'))