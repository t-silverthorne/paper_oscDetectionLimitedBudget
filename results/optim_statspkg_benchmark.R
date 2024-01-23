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
ampval

freqs=seq(from=1,to=24,length.out=1e3)
cfun = function(mt){1-costfun_svdpower(mt,freqs,Amp=2)}

cfun(runif(10))

# pwoer optimization on lattice
lat1=c(0:5)/5
lat2=c(0:5)/5
cfun = function(xpar){costfun_2lattice_svdpower(xpar[1],xpar[2],xpar[3],xpar[4],lat1,lat2,
                                         freqs,Amp=2)}
cfun(runif(4))