require(dplyr)
require(gurobi)
source('results/source/powerChord.R')

out_loc_glob = 'results/data/multi_highres/' 
Nmeas_vals   = seq(16,48,2)
w_reg_vals   = c(1,0)
drts_vals    = c(Inf,1*6,2*6,3*6)

setgs        = expand.grid(Nmeas=Nmeas_vals,w_reg=w_reg_vals,drts=drts_vals)
ii           = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
num_threads  = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

sol=powerChord(Nmeas       = setgs[ii,]$Nmeas,
               drts        = setgs[ii,]$drts,
               w_reg       = setgs[ii,]$w_reg,
               Nfreq       = 256,
               num_threads = num_threads,
               tlim        = 60*60)

saveRDS(sol,file=paste0(out_loc_glob,'sol_gur_',
                        'wreg_',setgs[ii,]$w_reg,
                        '_drts_',setgs[ii,]$drts,
                        '_Nmeas_',setgs[ii,]$Nmeas,'.RDS'))
