# setup 
require(devtools)
require(ggplot2)
require(annmatrix)
require(ggplotify)
require(patchwork)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(latex2exp)
source('plotConfig.R')
load_all()
fs_glob=9
sz_glob=1
mc_cores=10 # TODO: read this in from slurm

Nmeas = 16
Nmc   = 1e5 # TODO: read this in from slurm
Nfreq = Nmc
Nacro = Nmc
Amp   = 3
fmin  = 1
fmax  = 8
acros = 2*pi*runif(Nacro)
freqs = runif(Nfreq,fmin,fmax)
pars  = data.frame(Amp=Amp,acro=acros,freq=freqs)
tvec  = c(1:Nmeas)/Nmeas-1/Nmeas
monteCarloPval <-function(tvec,Amp,freq,acro){
  Ydat = matrix(Amp*cos(2*pi*freq*tvec -acro) + rnorm(length(tvec)),nrow=1)
  return(rowCosinor(Ydat,tvec,per=1/freq) %>% {.$pvalue} )
}

# Panel 1 of fig
pars$pval = c(1:dim(pars)[1]) %>%mclapply(mc.cores=mc_cores,function(ii){
  x=pars[ii,] 
  return(monteCarloPval(tvec,x[['Amp']],x[['freq']],x[['acro']]))
}) %>% unlist()
cticks=c(1,2,4,6,8)
pars = pars %>% mutate(Q=cut(freq,quantile(freq,0:4/4),include.lowest=T,include.heighest=T))
p1=pars %>% ggplot(aes(x=acro,y=pval,color=freq))+geom_point(size=sz_glob)+
  scale_color_viridis_c(limits=c(fmin,fmax),breaks=cticks,labels=24/cticks)+
  ylim(c(0,1))+
  scale_x_continuous(limits =c(0,2*pi),
                     breaks=c(0,pi/2,pi,3*pi/2,2*pi),
                     labels=c(TeX('$0$'),TeX('$\\pi/2$'),TeX('$\\pi$'),TeX('$3\\pi/2$'),TeX('$2\\pi$')))
saveRDS(p1,'figs_for_paper/Fig0_motivation/p1.RDS')

# Panel 2 of fig
NfreqFDR=2^5
freqsFDR = seq(from=fmin,to=fmax*.99,length.out=NfreqFDR)
tvec    = c(1:Nmeas)/Nmeas-1/Nmeas
acros   = 2*pi*runif(Nacro)
freqs   = runif(Nfreq,fmin,fmax)
for (jj in c(1:2)){
  if (jj ==1){
    tvec    = c(1:Nmeas)/Nmeas-1/Nmeas
  }else if (jj==2){
    set.seed(111)
    tvec    = runif(Nmeas)
  }else{
    stop('unrec ind')
  }
  pars=data.frame(Amp=Amp,acro=acros,freq=freqs)
  tmat    = replicate(Nmc,{tvec}) %>% t()
  acromat = replicate(Nmeas,{acros}) 
  X       = Amp*cos(2*pi*freqs*tmat-acromat)+matrix(rnorm(Nmeas*Nmc),nrow=Nmc)
  
  Lfit = freqsFDR %>% mclapply(mc.cores=mc_cores,function(freq){
    rowCosinor(X,zts=tvec,per=1/freq) %>% {data.frame(pvalue=t(.$pvalue))}
    })
      
  pmat = Lfit %>% rbindlist() %>% as.matrix() %>% t()  
  qvec = NULL
  for (ii in c(1:Nmc)){
        qvec[ii]=p.adjust(pmat[ii,],'BH') %>% min()
  }
  pars$qval = qvec
  
  if (jj==1){
    fdr_unif = pars 
  }else if(jj==2){
    fdr_rand = pars
  }else{
    stop('unrec ind')
  }
}

fdr_unif$sched='uniform'
fdr_rand$sched='random'
df=rbind(fdr_unif,fdr_rand)
df$sched=factor(df$sched,c('uniform','random'),c('uniform','random'))
df = df %>% mutate(Q=cut(freq,quantile(freq,0:4/4),include.lowest=T,include.heighest=T))

p2= df %>% ggplot(aes(x=acro,y=qval,color=freq))+geom_point(size=sz_glob)+
  scale_color_viridis_c(limits=c(fmin,fmax))+facet_grid(sched~Q)+
  scale_x_continuous(limits =c(0,2*pi),
                     breaks=c(0,pi/2,pi,3*pi/2,2*pi),
                     labels=c(TeX('$0$'),TeX('$\\pi/2$'),TeX('$\\pi$'),TeX('$3\\pi/2$'),TeX('$2\\pi$')))
saveRDS(p2,'figs_for_paper/Fig0_motivation/p2.RDS')