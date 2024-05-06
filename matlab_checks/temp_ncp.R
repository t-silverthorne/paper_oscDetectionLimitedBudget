require(parallel)
f = seq(1,4,.01)
t = seq(0,3,.01)

phi=1
pars=expand.grid(f=f,t=t)
h=1e-8
pars$fun_FD =c(1:dim(pars)[1]) %>% mclapply(function(ind){
  x=pars[ind,]
  f=as.numeric(x$f)
  t=as.numeric(x$t)
  funh = cos(2*pi*f*(1+h)*t)^2
  fun  = cos(2*pi*f*t)^2
  return((funh-fun)/f/h)
})

pars$fun_exact =c(1:dim(pars)[1]) %>% mclapply(function(ind){
  x=pars[ind,]
  f=as.numeric(x$f)
  t=as.numeric(x$t)
  return(4*t*pi*cos(phi - 2*pi*f*t)*sin(phi - 2*pi*f*t))
})
pars[1008,]

pars$fun_FD=as.numeric(pars$fun_FD)
pars$fun_exact=as.numeric(pars$fun_exact)
pars %>% ggplot(aes(x=t,y=f,fill=fun_exact))+geom_raster()+scale_fill_viridis_c()
pars %>% ggplot(aes(x=t,y=f,fill=fun_FD))+geom_raster()+scale_fill_viridis_c()
