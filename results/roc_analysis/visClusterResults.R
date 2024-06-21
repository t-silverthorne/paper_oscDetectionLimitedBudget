require(ggplot2)
require(dplyr)
require(patchwork)
df=readRDS('results/roc_analysis/worst_case_auc.RDS')
df$type = as.character(df$type)
df[df$type=='threshold',]$type    = 'WCP'
df[df$type=='balanced',]$type     = 'RVWCP constrained'
df[df$type=='regu_no_cstr',]$type = 'RVWCP free'
df$type=factor(df$type,levels=c('equispaced','WCP','RVWCP free','RVWCP constrained'))
dfw=df

df=readRDS('results/roc_analysis/revisedROC.RDS')
df$type = as.character(df$type)
df[df$type=='threshold',]$type    = 'WCP'
df[df$type=='balanced',]$type     = 'RVWCP constrained'
df[df$type=='regu_no_cstr',]$type = 'RVWCP free'
df$type=factor(df$type,levels=c('equispaced','WCP','RVWCP free','RVWCP constrained'))
head(df)
p1=df %>%mutate(N=Nmeas,A=Amp) %>%   ggplot(aes(x=freq,y=AUC,color=type,group=type))+
  geom_line()+
  facet_grid(A~N,labeller=purrr::partial(label_both,sep=' = '))


p2=dfw %>%mutate(N=Nmeas,A=Amp) %>%   ggplot(aes(x=freq,y=wAUC,color=type,group=type))+
  geom_line()+
  facet_grid(A~N,labeller=purrr::partial(label_both,sep=' = '))

p1 = p1 + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
p2 = p2 + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
p1/p2+ plot_annotation(tag_levels= 'A')+plot_layout(guides = 'collect')
