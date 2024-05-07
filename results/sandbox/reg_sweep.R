require(dplyr)
require(gurobi)
require(ggplot2)
require(data.table)
require(parallel)
source('results/source/powerChord.R')

tlim_glob  = 5 
fmin       = 1
fmax_vals  = seq(12,50,4)
#scale_vals = seq(1,1.5,.5)
Nmeas_vals = seq(40,60,5)
cfactor    = c(1.2,1.5)
pars = expand.grid(fmax=fmax_vals,Nmeas=Nmeas_vals,cfactor=cfactor)

pdf = c(1:dim(pars)[1]) %>% mclapply(mc.cores=10,function(ii){
  # unpack
  x       = pars[ii,]
  fmax    = x$fmax
  cfactor = x$cfactor 
  Nmeas   = x$Nmeas 
  Nfine   = floor(Nmeas*cfactor)
  
  sol=powerChord(Nmeas       = Nmeas,
                 drts        = Inf,
                 w_reg       = 1,
                 Nfreq       = NaN,
                 num_threads = 1,
                 tlim        = tlim_glob,
                 Nfine       = Nfine,
                 w_wc        = 0,
                 fmin        = fmin,
                 fmax        = fmax,
                 MIPGap      = 1e-12
  )   
    
  is_unif = all(as.logical(sol$x[1:Nfine]>.99) == 
      as.logical(c(rep(1,Nmeas),rep(0,Nfine-Nmeas))))
  return(data.frame(fmax=fmax,Nmeas=Nmeas,Nfine=Nfine,cfactor=cfactor,
                    status=sol$status,is_unif=is_unif))
  }
) %>% rbindlist()
saveRDS(pdf,'results/data/pure_reg.RDS')
#tpdf = pdf
#tpdf = tpdf %>% mutate(outcome = ifelse(status=='TIME_LIMIT','time_out',is_unif))
#tpdf$outcome = factor(tpdf$outcome,c(T,F,'time_out'),c('uniform optimal','non-uniform optimal','time limit'))
#tpdf %>% ggplot(aes(x=Nmeas,y=fmax,fill=outcome))  + geom_tile()+facet_wrap(~cfactor)
