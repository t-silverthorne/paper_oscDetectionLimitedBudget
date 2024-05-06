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



# for scale 1, the cutoff happens before fmax=12
# for scale 1.15, the cutoff happens after fmax=12
# for scale 1.45, the cutoff happens after fmax=14
scale=1.45
Nfine = floor(144*scale/3)
Nmeas = floor(40*scale)
fmin  = 1
fmax  = 14
sol=powerChord(Nmeas       = Nmeas,
               drts        = Inf,
               w_reg       = 1,
               Nfreq       = 8,
               num_threads = 12,
               tlim        = 35,
               Nfine       = Nfine,
               w_wc        = 0,
               fmin        = fmin,
               fmax        = fmax,
               MIPGap      = 1e-12
               ) 

tau   = c(1:Nfine)/Nfine -1/Nfine 

data.frame(time=tau[sol$x[1:Nfine]>0]) %>% 
  ggplot(aes(x=time,y=1))+geom_point()+
  scale_x_continuous(limits=c(0,1))

#compare costfun
mt_opt = tau[sol$x[1:Nfine]>0]
length(mt_opt)==Nmeas

mt_tight= tau[1:Nmeas]

mt_unif = c(1:Nmeas)/Nmeas-1/Nmeas

get_exact_reg = function(tau,fmin,fmax){
  dtau_mat          = outer(tau,tau,'-')
  tau_prod_mat      = outer(tau,tau,'*')
  ncp_reg_mat       = ((sin(4*pi*fmax*dtau_mat) - sin(4*pi*fmin*dtau_mat))*(4*tau_prod_mat*pi^2 + 1))/(4*dtau_mat)
  diag(ncp_reg_mat) = pi*(4*pi^2*tau^2 + 1)*(fmax - fmin) 
  return(sum(sum(ncp_reg_mat)))
}

mt_rand = runif(Nmeas)
get_exact_reg(mt_rand,1,24)
get_exact_reg(mt_rand+runif(1),1,24)


get_exact_reg(mt_opt,fmin,fmax)
get_exact_reg(mt_tight,fmin,fmax)
get_exact_reg(mt_unif,fmin,fmax)
get_exact_reg(tau[sample(c(1:Nfine),Nmeas)],fmin,fmax)

