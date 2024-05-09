require(lomb)
Nmeas = 48
Nmc   = 1e3
Amp   = 1 
freq  = 5.05

estimate_TPR = function(alpha_val,p_osc,tvec,Nmc,Amp,freq,type){
  # generate simulated data
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
 
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb           = lsp(x,alpha=alpha_val,tvec,plot=F)
    freq_est       = lomb$peak.at[1]
    pcosinor       = rowCosinor(matrix(x,nrow=1),matrix(tvec,ncol=1),freq_est) %>% {.$pvalue}
    return(data.frame(plomb=lomb$p.value,pcosinor=pcosinor))
  }) %>% rbindlist() %>% data.frame()
  num_P  = sum(state=='osc')
  
  num_TP1 = sum(pvdf$pcosinor<alpha_val & state=='osc')
  num_TP2 = sum(pvdf$plomb <alpha_val & state=='osc')
  num_TP3 = sum(pvdf$pcosinor<alpha_val &pvdf$plomb <alpha_val & state=='osc')
  return(data.frame(alpha_val=alpha_val,
                    freq=freq,
                    type=type,
                    TPR1=num_TP1/num_P,
                    TPR2=num_TP2/num_P,
                    TPR3=num_TP3/num_P))
}
alpha_val = 0.5 
p_osc     = 0.1 


alpha_vals = seq(0,1,.1)
freq_vals  = c(1,2,4,24) 
mt_unif    = c(1:Nmeas)/Nmeas-1/Nmeas
mt_opt     = sols[sols@wreg==0 & sols@drts==Inf & sols@Nmeas==Nmeas,]
pars       = expand.grid(alpha_val=alpha_vals,freq=freq_vals,type=c('uniform','optimal'))

dfTPR = c(1:dim(pars)[1]) %>% lapply(function(ind){
  alpha_val = pars[ind,]$alpha_val
  freq      = pars[ind,]$freq
  if (pars[ind,]$type=='uniform'){
    tvec_loc = mt_unif
  }else{
    tvec_loc = mt_opt
  }
  estimate_TPR(alpha_val,p_osc,tvec_loc,Nmc,Amp,freq,pars[ind,]$type)
}
) %>%  rbindlist() %>% data.frame()
pdf =dfTPR %>% gather(key='TPR_type',value='TPR',TPR1,TPR2,TPR3)
pdf %>% ggplot(aes(x=alpha_val,y=TPR,group=type,color=type))+geom_line()+
  scale_y_continuous(limits = c(0,1))+
  facet_grid(TPR_type~freq)