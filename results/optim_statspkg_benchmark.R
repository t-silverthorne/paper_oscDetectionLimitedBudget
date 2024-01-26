library(devtools)
load_all('.')
library(stats)
cfun = function(x){sum((x-1.23)^2)}
d=10
start=runif(d)
lower=rep(0,d)
upper=rep(10,d)
xopt=optim(start,cfun,gr=function(x){x+rnorm(length(x))},
              method='SANN')
xopt

# power optimization
control = list()
control$Nmeas       = 16
control$ampval      = 2
control$tfun_choice = 'brownian'
control$mean        = 0
control$sd          = .01
freqs = seq(from=1,to=24,length.out=1e3)
cfun  = function(mt){1-costfun_svdpower(mt,freqs,Amp=control$ampval)}
tfun  = function(mt){tfun_svdpower(mt,control)}

mt0 = runif(control$Nmeas)
xopt=stats::optim(mt0,fn=cfun,gr=tfun,
              method='SANN',control=list(trace=100,REPORT=1,maxit=100))

xopt=stats::optim(mt0,fn=cfun,gr=NULL,
                  lower=rep(0,length(mt0)),
                  upper=rep(1,length(mt0)),
              method='L-BFGS-B',control=list(trace=6,REPORT=1,maxit=1))


# pwoer optimization on lattice
N1=10
N2=6
lat1=c(0:(N1-1))/N1
lat2=c(0:(N2-1))/N2
cfun = function(xpar){1-costfun_2lattice_svdpower(0,xpar[1],xpar[2],xpar[3],lat1,lat2,
                                         freqs,Amp=2)}
lat2_lower = c(0,1e-3,1e-3)
lat2_upper = c(1/N2,1,1)

start=Sys.time()
#mt0 = c(1,1,0,1/N) 
mt0 =runif(3) 
xopt=stats::optim(mt0,fn=cfun,gr=NULL,
                  #lower=lat2_lower,
                  #upper=lat2_upper,
                  method='L-BFGS-B',control=list(trace=6,REPORT=1,maxit=2))
end=Sys.time()
end-start

pts=convert_2lattice_to_state(0,xopt$par[1],xopt$par[2],xopt$par[3],lat1,lat2)
plot(pts)
