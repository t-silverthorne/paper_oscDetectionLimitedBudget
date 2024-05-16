require(lomb)
require(parallel)
require(data.table)
#' Estimate TPR for multi-frequency hypothesis testing
#' 
#' @param alpha the false positive rate
#' @param Nmc number of Monte-Carlo replicates
#' @param tvec vector of measurement times
#' @param Amp amplitude of signal
#' @param freq frequency of signal
#' @param p_osc what proportion of the simulated data contains oscillations
#' @param mc_cores number of cores to be used for parallelization
#' @return a dataframe containing estimate of the TPR for each combination FDR method
#' and hyptest given in \code{FDR_list} and \code{hyptest_list}
#' 
#' @author Turner Silverthorne
#' @export 
estimateTPR=function(alpha_vals,Nmc,tvec,
                       Amp=1,freq=1,p_osc=.1,mc_cores=1,
                     fdr_methods=c('none','BH')){
  # simulated data
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
  
  # test periodogram signif on each signal 
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,alpha=alpha,tvec,plot=F,normalize = 'standard')
    lomb_press     = lsp(x,alpha=alpha,tvec,plot=F,normalize = 'press')
    return(rbind(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]),
                 data.frame(p_method='press',pval =lomb_press$p.value,state=state[ii])))
  }) %>% rbindlist() %>% data.frame()
  num_P  = sum(state=='osc')
  num_N  = sum(state!='osc')
 
  if (all(fdr_methods  %in% c('none','BH','BY'))){
    methods = expand.grid(p_method = c('std','press'),fdr_method=fdr_methods)
  }else{
    stop('unknown fdr method')
  }
  
  methods$TPR = methods %>% apply(1,function(x){
    dfloc = pvdf[pvdf$p_method==x[['p_method']],]
    if (x[['fdr_method']]=='none'){
      ploc = dfloc$pval 
    }else{
      ploc = p.adjust(dfloc$pval,method=x[['fdr_method']]) 
    }
    return(sum(ploc<alpha & state == 'osc')/num_P)
  })
  
  methods$FPR = methods %>% apply(1,function(x){
    dfloc = pvdf[pvdf$p_method==x[['p_method']],]
    if (x[['fdr_method']]=='none'){
      ploc = dfloc$pval 
    }else{
      ploc = p.adjust(dfloc$pval,method=x[['fdr_method']]) 
    }
    return(sum(ploc<alpha & state == 'non_osc')/num_N)
  })
      
  return(methods)  
}