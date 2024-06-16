rm(list=ls(all=T)) # clear env
source('results/figure_generation/plos_figs/fig_env.R')
require(ggplot2)
require(patchwork)
require(dplyr)
require(data.table)
require(devtools)
require(purrr)
load_all()
loadedNamespaces()

fsize=9
theme_set(theme_classic())

cmap = c('RVWCP'     ='#40826E',
         'WCP'       ='#A60A5E',
         'random'    ='#003399',
         'equispaced'='black'
)

df=readRDS('results/data/roc.RDS')
df$type = as.character(df$type)
df[df$type=='threshold',]$type    = 'WCP'
df[df$type=='balanced',]$type     = 'RVWCP constrained'
df[df$type=='regu_no_cstr',]$type = 'RVWCP free'
df$type=factor(df$type,levels=c('equispaced','WCP','RVWCP free','RVWCP constrained'))
df.sum=df %>% filter(type!='random' ) %>%
  group_by(freq,Nmeas,Amp,p_osc,fdr_method,type) %>% 
  summarise(sd_AUC = sd(AUC),
            sd_TPR = sd(TPR),
            sd_FPR = sd(FPR),
            AUC = mean(AUC),
            TPR=mean(TPR),
            FPR=mean(FPR))

plt=df.sum %>% mutate(N=Nmeas,A=Amp) %>% 
  ggplot(aes(x=freq,y=AUC,color=type,group=type))+
  geom_line()+#geom_errorbar(aes(ymin=AUC-sd_AUC,ymax=AUC+sd_AUC),data=df.sum)+
  facet_grid(N~A,labeller = purrr::partial(label_both,sep=' = '))+
  scale_color_brewer(type='qual',palette =6)+
  scale_x_continuous(limits = c(1,24),breaks=c(1,seq(4,24,4)))+
  labs(x=element_text('frequency'))

plt=plt+ theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt_height = 4
plt_width  = 6
plt=plt+theme(legend.position='bottom',
               legend.box='vertical',
               legend.direction = "horizontal",
               legend.margin=margin())


plt=plt+theme(text=element_text(size=fsize))
show_temp_plt=function(plt,plt_width,plt_height){
  plt_path <- tempfile(fileext = ".png")
  ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
         dpi = 96)
  
  viewer <- getOption("viewer")
  viewer(plt_path)
}

show_temp_plt(plt,plt_width,plt_height)
ggsave(paste0('~/research/ms_powerCHORD/figures/plos_submission/','fig5.png'),plt,width=plt_width,height=plt_height,device='png',dpi=600)
