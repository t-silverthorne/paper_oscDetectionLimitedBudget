require(dplyr)
require(gurobi)
source('results/source/powerChord.R')

out_loc_glob = 'results/data/' 
Nmeas_vals=c(24:32)
w_reg_vals=c(1e-5,1e-2,1,0)
drts_vals =c(Inf,1*6,2*6,3*6)
setgs =expand.grid(Nmeas=Nmeas_vals,w_reg=w_reg_vals,drts=drts_vals)
ii=1
sol=powerChord(Nmeas=setgs[ii,]$Nmeas,
               drts=setgs[ii,]$drts,
               w_reg=setgs[ii,]$w_reg,
               num_threads=1,
               tlim=60*60)

saveRDS(sol,file=paste0(out_loc_glob,'sol_gur_',
                        'wreg_',w_reg,
                        '_drts_',drts,
                        '_Nmeas_',Nmeas,'.RDS'))
