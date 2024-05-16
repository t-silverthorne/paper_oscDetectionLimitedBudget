source('figs_for_paper/Fig1_aug_vs_cvxr/settings.R')
source('plotConfig.R')
mc_cores=Sys.getenv("SLURM_CPUS_PER_TASK") 

# load solutions
mtdf_full=readRDS('figs_for_paper/Fig1_aug_vs_cvxr/solns.RDS')
head(mtdf_full)

# make panel 1
p1=mtdf_full %>% ggplot(aes(x=time,y=ncp,color=solver))+xlim(c(0,1))+geom_point(size=.7)+
  facet_wrap(~Nmeas,ncol=length(Nmeasvals))
saveRDS(p1,'figs_for_paper/Fig1_aug_vs_cvxr/p1.RDS')

# make panel 2
amps = seq(from=.5,3.0,length.out=N_amp_plt)
pers = 10^seq(from=log10(1/24),to=log10(1),length.out=N_per_plt)
tag   = c('uniform','pin2lat','cvxr')
pars  = expand.grid(amp=amps,per=pers,tag=tag,Nmeas=Nmeasvals)
minPower = rep(NaN,dim(pars)[1])
minPower =c(1:length(minPower)) %>% as.list() %>% mclapply(mc.cores=mc_cores,
  function(ii){
    x=pars[ii,]
    if (x[['tag']]=='uniform'){
      Nmeas  = as.numeric(x[['Nmeas']])
      mt_loc = c(1:Nmeas)/Nmeas-1/Nmeas
    }else{
      filt   = mtdf_full$Nmeas==as.numeric(x[['Nmeas']]) & mtdf_full$solver==x[['tag']]
      mt_loc = as.numeric(mtdf_full[filt,]$time)
    }
      return(costfun_svdpower(mt    = mt_loc,
                              freqs = 1/as.numeric(x[['per']]),
                              Amp   = as.numeric(x[['amp']])))
  }
)%>% unlist()
pars$minPower=minPower
pars$circ_per = 24*pars$per
p2=pars %>% ggplot(aes(x=amp,y=circ_per,fill=minPower))+geom_tile(show.legend=F)+facet_grid(tag~Nmeas)+
  scale_fill_viridis_c(limits=c(0,1))+
  scale_x_continuous(breaks=seq(min(amps),max(amps),length.out=3))+
  scale_y_continuous(breaks=c(c(1:2),seq(4,24,4)),trans='log10')
p2=p2+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
saveRDS(p2,'figs_for_paper/Fig1_aug_vs_cvxr/p2.RDS')




# make panel 3
amps = seq(from=.5,3.0,length.out=N_amp_plt)
pers = 10^seq(from=log10(1/24),to=log10(1),length.out=N_per_plt)
tag   = c('uniform','pin2lat','cvxr')
fdrs  = expand.grid(amp=amps,per=pers,tag=tag,Nmeas=Nmeasvals)
minFDR = replicate(dim(fdrs)[1],{NaN})
minFDR= c(1:dim(fdrs)[1]) %>% as.list() %>% 
  mclapply(mc.cores=mc_cores,function(ii){
  x=fdrs[ii,]
  if (x[['tag']]=='uniform'){
    Nmeas  = as.numeric(x[['Nmeas']])
    mt_loc = c(1:Nmeas)/Nmeas-1/Nmeas
  }else{
    filt   = mtdf_full$Nmeas==as.numeric(x[['Nmeas']]) & mtdf_full$solver==x[['tag']]
    mt_loc = as.numeric(mtdf_full[filt,]$time)
  }
    return(benchmark_worst_fdr(mt = mt_loc,
                            f0    = 1/as.numeric(x[['per']]),
                            Amp   = as.numeric(x[['amp']]),
                            Nmc   = NmcFDR,
                            Nacro = NacroFDR,
                            Nfreq = NfreqFDR,
                            fmin  = fmin,
                            fmax  = fmax*.99)
           )
  }) %>% unlist()

fdrs$minFDR = as.numeric(minFDR)
fdrs$circ_per=24*fdrs$per
p3=fdrs %>% ggplot(aes(x=amp,y=circ_per,fill=minFDR))+geom_tile()+facet_grid(tag~Nmeas)+
  scale_fill_viridis_c(limits=c(0,1))+
  scale_x_continuous(breaks=seq(min(amps),max(amps),length.out=3))+
  scale_y_continuous(breaks=c(c(1:2),seq(4,24,4)),trans='log10')
  labs(fill = 'power',
          y = 'frequency',
          x = 'amplitude')
p3=p3+theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())
saveRDS(p3,'figs_for_paper/Fig1_aug_vs_cvxr/p3.RDS')