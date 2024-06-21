require(lomb)
Nmeas = 48
Nmc   = 1e4
Amp   = 1 

estimate_TPR = function(alpha_val,p_osc,tvec,Nmc,Amp,freq,type){
  # generate simulated data
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
 
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,alpha=alpha_val,tvec,plot=F,normalize = 'standard')
    lomb_press     = lsp(x,alpha=alpha_val,tvec,plot=F,normalize = 'press')
    freq_est       = lomb$peak.at[1]
    pcosinor       = row_cosinor(matrix(x,nrow=1),matrix(tvec,ncol=1),freq_est) %>% {.$pvalue}
    return(data.frame(plomb_std=lomb_std$p.value,
                      plomb_press=lomb_press$p.value,
                      pcosinor=pcosinor))
  }) %>% rbindlist() %>% data.frame()
  num_P  = sum(state=='osc')
  
  num_TP1 = sum(pvdf$pcosinor<alpha_val & state=='osc')
  num_TP2 = sum(pvdf$plomb_std <alpha_val & state=='osc')
  num_TP3 = sum(pvdf$plomb_press <alpha_val & state=='osc')
  return(data.frame(alpha_val=alpha_val,
                    freq=freq,
                    Amp =Amp,
                    type=type,
                    TPR1=num_TP1/num_P,
                    TPR2=num_TP2/num_P,
                    TPR3=num_TP3/num_P))
}
alpha_val = 0.5 
p_osc     = 0.1 


alpha_vals = seq(0,1,.1)
freq_vals  = c(1,2,4,8,16,24) 
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
pars       = expand.grid(alpha_val=alpha_vals,
                         freq=freq_vals,
                         Amp = c(1,2,10),
                         type=c('uniform','optimal','robust'))

dfTPR = c(1:dim(pars)[1]) %>% lapply(function(ind){
  alpha_val = pars[ind,]$alpha_val
  freq      = pars[ind,]$freq
  Amp       = pars[ind,]$Amp
  if (pars[ind,]$type=='uniform'){
    tvec_loc = mt_unif
  }else if(pars[ind,]$type=='optimal'){
    tvec_loc = mt_opt
  }else if(pars[ind,]$type=='robust'){
    tvec_loc = mt_rob
  }else{
    stop('unrecognized type')
  }
  estimate_TPR(alpha_val,p_osc,tvec_loc,Nmc,Amp,freq,pars[ind,]$type)
}
) %>%  rbindlist() %>% data.frame()

pdf =dfTPR %>% gather(key='TPR_type',value='TPR',TPR1,TPR2,TPR3)
saveRDS(pdf,'results/data/hilo_ROC.RDS')
pdf = readRDS('results/data/hilo_ROC.RDS')
plt1 =pdf %>% filter(Amp==1)%>% ggplot(aes(x=alpha_val,y=TPR,group=type,color=type))+
  geom_line()+
  #geom_line(position=position_jitter(w=0.02, h=0.02))+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(TPR_type~freq)

plt2 =pdf %>% filter(Amp==2) %>% ggplot(aes(x=alpha_val,y=TPR,group=type,color=type))+
  geom_line()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(TPR_type~freq)

plt3 =pdf %>% filter(Amp==10) %>% ggplot(aes(x=alpha_val,y=TPR,group=type,color=type))+
  geom_line()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(TPR_type~freq)

theme_set(theme_classic()) 
plt_width  = 6 
plt_height = 7
plt = plt1/plt2

modify_plt=function(plt){
  plt = plt + labs(x=element_text('FPR'),
                   y=element_text('TPR'),
                   color='schedule')
  
  xmin=0
  xmax=1
  plt = plt + scale_x_continuous(limits =c(xmin,xmax),
                                 breaks=seq(xmin,xmax,.25),
                                 labels=seq(xmin,xmax,.25)) 
  plt = plt + scale_y_continuous(limits =c(xmin,xmax),
                                 breaks=seq(xmin,xmax,.25),
                                 labels=seq(xmin,xmax,.25)) 
  
  #plt=plt+guides(color=guide_colorbar(title.position='top'))
  #plt=plt+theme(legend.position='bottom',
  #              legend.key.width = unit(plt_width*.05, "in"),
  #              legend.title.align = 0.5,
  #              legend.direction = "horizontal")
  plt=plt+theme(text=element_text(size=9))
  return(plt)
}
plt1_mod = modify_plt(plt1)
plt2_mod = modify_plt(plt2)


Fig = (plt1_mod/plt2_mod)+plot_annotation(tag_levels='A') + plot_layout(guides='collect') &
  theme(legend.position='bottom',
               legend.key.width = unit(plt_width*.05, "in"),
               legend.title.align = 0.5,
               legend.direction = "horizontal")
Fig = Fig & theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

show_temp_plt(Fig,plt_width,plt_height)
ggsave(paste0('~/research/ms_powerCHORD/figures/','roc_analysis.png'),Fig,
       width=plt_width,height=plt_height,device='png',dpi=600)