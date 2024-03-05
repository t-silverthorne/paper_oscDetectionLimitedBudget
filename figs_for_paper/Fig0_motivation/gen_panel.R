gen_panel=function(mt){
  Nmeas =length(mt)
  mtdf=data.frame(time=mt)
  
  # panel 1 
  t_fine = seq(0,1,.005)
  freqs = c(1,6,16)
  cp_levels = sort(24/freqs) %>% rev()
  cp_labels = circ_per_levels %>%  {paste0(.,' hr')}
  phi_plt1=pi/2
  
  df = freqs %>% lapply(function(freq){
   data.frame(expr=cos(2*pi*freq*t_fine - phi_plt1),freq=freq,time=t_fine)
  }) %>% rbindlist() %>% data.frame()
  
  df_noise = freqs %>% lapply(function(freq){
   data.frame(expr=cos(2*pi*freq*mt - phi_plt1)+rnorm(length(mt)),
              freq=freq,time=mt)
  }) %>% rbindlist() %>% data.frame()
  df$circ_per=24/df$freq
  df$circ_per = factor(df$circ_per,cp_levels,cp_labels)
  df_noise$circ_per=24/df_noise$freq
  df_noise$circ_per = factor(df_noise$circ_per,cp_levels,cp_labels)
  p1=ggplot()+geom_line(data=df,aes(x=24*time,y=expr,group=circ_per,color=circ_per))+
    scale_color_manual(values = freq_colors)+
    geom_point(data=df_noise,aes(x=24*time,y=expr,color=circ_per),size=.5)+facet_wrap(~circ_per,nrow=3)+labs(y='expression',x='time (hr)')+
    scale_x_continuous(limits=c(0,24),breaks = 4*c(0:6))
 
  # panel 2 
  Nmc=1e2
  Amp=2
  monteCarloPval <-function(tvec,Amp,freq,acro){
    Ydat = replicate(Nmc,{Amp*cos(2*pi*freq*tvec -acro) + rnorm(length(tvec))}) %>% t
    return(rowCosinor(Ydat,tvec,per=1/freq) %>% {.$pvalue} )
  }
  acrovec=seq(from=0,to=2*pi,.05)
  pars=expand.grid(freq=freqs,acro=acrovec)
  
  #TODO: better way of doing this
  df=c(1:dim(pars)[1]) %>% lapply(function(ind){
    x=pars[ind,]
    freq=as.numeric(x['freq'])
    acro=as.numeric(x['acro'])
    return(data.frame(freq=freq,acro=acro,
                      pvalue=monteCarloPval(mt,Amp,freq*.999,acro)))
  }) %>% rbindlist() %>% data.frame()
  
  
  df$circ_per=24/df$freq
  df$circ_per = factor(df$circ_per,cp_levels,cp_labels)
  
  p2=df %>% ggplot(aes(x=acro,y=pvalue,color=circ_per))+geom_point(size=.5,alpha=.2)+facet_wrap(~circ_per,nrow=3)+
    scale_y_continuous(trans='log10')+
    geom_hline(aes(yintercept =.05,linetype='p=0.05'),color='black')+
    scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk,labels = rad_lab)+labs(x='acrophase (rad)',linetype=element_blank())+
    scale_color_manual(values=freq_colors)+scale_linetype_manual(values=c('dashed'))+guides(color='none')

  # panel 3
  acros = seq(0,2*pi,.01)
  pars=expand.grid(freq=freqs,acro=acros) 
  df= c(1:dim(pars)[1]) %>% lapply(function(ind){
    x=pars[ind,]
    param=list(
      Amp=Amp,
      freq=as.numeric(x['freq']),
      acro=as.numeric(x['acro'])
    )
    return(data.frame(freq=param$freq,acro=param$acro,power=eval_exact_power(mt,param)))
  }) %>% rbindlist() %>% data.frame()
  
  df$circ_per = 24/df$freq
  df$circ_per = factor(df$circ_per,cp_levels,cp_labels)
  
  p3=df %>% ggplot(aes(x=acro,y=power,color=circ_per))+geom_line()+facet_wrap(~circ_per,nrow=3)+labs(x='acrophase (rad)')+
    scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk,labels = rad_lab)+
    scale_color_manual(values=freq_colors)+guides(color='none')+ylim(c(0,1))

  Fig=p1+p2+p3  + plot_layout(guides='collect') & theme(
     strip.background=element_blank(),
     text=element_text(size=9),
     plot.margin=margin(0,0,0,0),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     axis.text.x = element_text(vjust = 0.5)
  )& labs(color='period')
 return(Fig) 
}