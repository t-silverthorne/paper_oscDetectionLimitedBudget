df=readRDS('results/data/roc2_hires.RDS')

df.sum=df %>% filter(type!='random' & fdr_method=='fdr') %>% group_by(Amp,p_osc,type,freq) %>% 
  summarise(sd_AUC = sd(AUC),
            sd_TPR = sd(TPR),
            sd_FPR = sd(FPR),
            AUC = mean(AUC),
            TPR=mean(TPR),
            FPR=mean(FPR))
df.sum %>%  ggplot(aes(x=freq,y=AUC,color=type,group=type))+
  geom_line()+geom_errorbar(aes(ymin=AUC-sd_AUC,ymax=AUC+sd_AUC),data=df.sum)+
  facet_grid(Amp~p_osc)#+ylim(c(0.25,1))


#df.sum %>% ggplot(aes(x=freq,y=TPR,color=type))+geom_line()+facet_grid(Amp~p_osc)
#df.sum %>% ggplot(aes(x=freq,y=FPR,color=type))+geom_line()+facet_grid(Amp~p_osc)
#
#df %>% group_by(Amp,p_osc,type) %>% summarise(amean=mean(AUC))
#
#for (ii in c(1:100)){
#  x              = Ydat[ii,]
#  lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
#}
