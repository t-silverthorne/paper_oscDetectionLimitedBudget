require(dplyr)
require(gurobi)
require(ggplot2)
source('results/source/powerChord.R')



## Example of non-uniform optimal solution
#Nfine = 144/4
#sol=powerChord(Nmeas       = 24,
#               drts        = Inf,
#               w_reg       = 1,
#               Nfreq       = Inf,
#               num_threads = 12,
#               tlim        = 20,
#               Nfine       = Nfine,
#               w_wc        = 0,
#               fmin        = 1,
#               fmax        = 7,
#               MIPGap      = 1e-12
#               ) 
#tau   = c(1:Nfine)/Nfine -1/Nfine 
#
#data.frame(time=tau[sol$x[1:Nfine]>0]) %>% 
#  ggplot(aes(x=time,y=1))+geom_point()+
#  scale_x_continuous(limits=c(0,1))

Nfine = 144 
sol=powerChord(Nmeas       = 24,
               drts        = Inf,
               w_reg       = 1,
               Nfreq       = Inf,
               num_threads = 12,
               tlim        = 20,
               Nfine       = Nfine,
               w_wc        = 0,
               fmin        = 12,
               fmax        = 24,
               MIPGap      = 1e-12
               ) 
tau   = c(1:Nfine)/Nfine -1/Nfine 

data.frame(time=tau[sol$x[1:Nfine]>0]) %>% 
  ggplot(aes(x=time,y=1))+geom_point()+
  scale_x_continuous(limits=c(0,1))