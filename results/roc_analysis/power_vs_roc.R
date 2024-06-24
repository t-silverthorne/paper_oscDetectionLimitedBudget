require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(lomb)
require(pROC)
require(ggplot2)
load_all()
Amp=1

sols = readRDS('results/data/MCperiodogram/hiresSols.RDS')

# Assemble AUC dataframe
df1        = readRDS('results/roc_analysis/revisedROC.RDS')
df2        = readRDS('results/roc_analysis/revisedROC_highfreq.RDS')
dfAUC      = rbind(df1,df2)
dfAUC      = dfAUC[dfAUC$Amp==Amp,]
colnames(dfAUC)[colnames(dfAUC)=='AUC']='FOM'
dfAUC$stat = 'Lomb-mean-AUC'
dfAUC      = dfAUC %>% select(!c(TPR,FPR,fdr_method,p_osc)) 

# Construct power dataframe
freqs = seq(1,36,.05)

pars       = expand.grid(freq=freqs,
                         Nmeas=c(32,40,48),
                         Amp = c(Amp),
                         type=c('equispaced','threshold','balanced','regu_no_cstr'))

dfpwr = c(1:dim(pars)[1]) %>% mclapply(mc.cores=12,function(ind){
  freq  = pars[ind,]$freq
  Nmeas = pars[ind,]$Nmeas
  
  mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
  mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
  mt_rob     = sols[sols@wreg==1 & sols@drts==6 & sols@Nmeas==Nmeas,]
  mt_rnc     = sols[sols@wreg==1 & sols@drts==Inf & sols@Nmeas==Nmeas,]
  mt_opt     = mt_opt[mt_opt<Inf]
  mt_rob     = mt_rob[mt_rob<Inf]
  mt_rnc     = mt_rnc[mt_rnc<Inf]
  
  if (pars[ind,]$type=='equispaced'){
    tvec = mt_unif
  }else if(pars[ind,]$type=='threshold'){
    tvec = mt_opt
  }else if(pars[ind,]$type=='balanced'){
    tvec = mt_rob
  }else if(pars[ind,]$type=='regu_no_cstr'){
    tvec = mt_rnc
  }else{
    stop('unrecognized type')
  }
  return(cbind(pars[ind,],data.frame(FOM=evalWorstPower(tvec,freq,Amp))))
}) %>% rbindlist() %>% data.frame()
dfpwr$stat='cosinor-min-power'

df = rbind(dfAUC,dfpwr)

# Get power AUC df
dfpAUC = readRDS('results/roc_analysis/cosinorROC.RDS')
dfpAUC = dfpAUC %>% select(!c(TPR,FPR,fdr_method,p_osc,stat_method)) 
colnames(dfpAUC)[colnames(dfpAUC)=='AUC']='FOM'
dfpAUC$stat='cosinor-mean-AUC'
df=rbind(df,dfpAUC)

# Get empirical TPR
dfcTPR= readRDS('results/roc_analysis/cosinorROC.RDS')
dfcTPR = dfcTPR %>% select(!c(FPR,AUC,fdr_method,p_osc,stat_method)) 
colnames(dfcTPR)[colnames(dfcTPR)=='TPR']='FOM'
dfcTPR$stat='cosinor-mean-TPR'

df=rbind(df,dfcTPR)

df$type = as.character(df$type)
df[df$type=='threshold',]$type    = 'WCP'
df[df$type=='balanced',]$type     = 'RVWCP constrained'
df[df$type=='regu_no_cstr',]$type = 'RVWCP free'
df$type=factor(df$type,levels=c('equispaced','WCP','RVWCP free','RVWCP constrained'))

df$stat = factor(df$stat,levels=c('cosinor-min-power','cosinor-mean-TPR','cosinor-mean-AUC','Lomb-mean-AUC'))
df %>% ggplot(aes(x=freq,y=FOM,group=type,color=type))+geom_line()+
  facet_grid(Nmeas~stat)

# Merge

# Make plot