require(data.table)
require(ggplot2)
require(dplyr)
require(parallel)

fmin=1
fmax=24
Nfreq=2^4
freqs = seq(fmin,fmax,length.out=Nfreq)
eps = 2*pi/Nmeas
acros = eps+seq(0,2*pi,length.out=2^2) 

lsp_pseudowindow=function(mt,freqs,acros,type,Nmeas){
  ploc = expand.grid(freq=freqs,acro=acros)
  c(1:dim(ploc)[1]) %>% lapply(function(ind){
    freq=ploc[ind,]$freq
    acro=ploc[ind,]$acro
    pl=lsp(x=cos(2*pi*freq*mt-acro),times=mt,from=fmin,to=fmax,plot=F)
    data.frame(lfreq=pl$scanned,lpower=pl$power,freq=freq,acro=acro,type=type,Nmeas=Nmeas)
  }) %>% rbindlist()
}

am=readRDS(file = 'figs_for_paper/Fig1_aug_vs_cvxr/solns.RDS')
pars=expand.grid(Nmeas=unique(am$Nmeas),tag=c('uniform',unique(am$solver)))

df=c(1:dim(pars)[1]) %>% mclapply(function(ii){
   x=pars[ii,]
    if (x[['tag']]=='uniform'){
      Nmeas  = as.numeric(x[['Nmeas']])
      mt_loc = c(1:Nmeas)/Nmeas-1/Nmeas
    }else{
      filt   = am$Nmeas==as.numeric(x[['Nmeas']]) & am$solver==x[['tag']]
      mt_loc = as.numeric(am[filt,]$time)
    }
  return(lsp_pseudowindow(mt_loc,freqs,acros,x[['tag']],as.numeric(x[['Nmeas']])))
}) %>% rbindlist()

head(df)

pdf = df[df$Nmeas==16,]
pdf %>% 
  ggplot(aes(x=freq,y=lfreq,fill=lpower))+geom_raster() +facet_grid(type~acro)+
  scale_fill_viridis_c()

