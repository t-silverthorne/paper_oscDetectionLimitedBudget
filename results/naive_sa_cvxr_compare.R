gc()
library(devtools)
load_all('.')
opts=make_default_opts(prob_size='medium',solver_type='simulanneal')

opts$Nfine = 1e3
opts$Nmeas = 16
opts$fmin  = 1
opts$fmax  = 24
opts$Nfreq = 2^9
opts$costfun_type = 'Linfty'
opts$verbose=T
opts$num_iter=150

rm(Aquad)
Aquad=make_quadmats(opts)

statime=Sys.time()
xopt=run_sa_power(Aquad,opts)
endtime=Sys.time()
endtime-statime

tau   = (0:opts$Nfine)/opts$Nfine       
tau   = tau[1:(length(tau)-1)] 

mt_opt    = tau[xopt$xval==1]

Nacro=2^8
Amp=1.5
freqs = seq(from=opts$fmin,to=opts$fmax,length.out=opts$Nfreq)
acros = seq(from=0,to=2*pi,length.out=Nacro+1)
acros = acros[1:Nacro]
pars  = expand.grid(Amp=Amp,freq=freqs,acro=acros)

worst_power = function(mt){
  min(apply(pars,1,function(x){eval_exact_power(t=mt,param=x)}))
}

#worst_power(c(0:10)/10)
1-worst_power(mt_opt)
