require(lomb)
require(parallel)
require(data.table)
require(dplyr)
#' Estimate TPR for multi-frequency hypothesis testing
#' 
#' @param alpha_vals values of threshold parameter
#' @param Nmc number of signals to include in dataset
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
getLombScargleROC=function(alpha_vals,Nmc,tvec,
                       Amp=1,freq=1,p_osc=.1,mc_cores=1,
                     fdr_methods=c('none','BH')){
  # generate simulated data
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
  num_P  = sum(state=='osc')
  num_N  = sum(state=='non_osc')
  
  # compute p-values from Lomb-Scargle
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
    lomb_press     = lsp(x,times=tvec,plot=F,normalize = 'press')
    return(rbind(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]),
                 data.frame(p_method='press',pval =lomb_press$p.value,state=state[ii])))
  }) %>% rbindlist() %>% data.frame()
 
  # check multiple test correction methods are implemented 
  if (all(fdr_methods  %in% c('none','BH','BY'))){
    methods = expand.grid(p_method = c('std','press'),fdr_method=fdr_methods,alpha=alpha_vals)
  }else{
    stop('unknown fdr method')
  }

  ROCdf =c(1:dim(methods)[1]) %>% lapply(function(ind){
    p_method   = toString(methods[ind,]$p_method)
    fdr_method = toString(methods[ind,]$fdr_method)
    alpha      = as.numeric(methods[ind,]$alpha)
    dfloc = pvdf[pvdf$p_method==p_method,]
    if (fdr_method=='none'){
      ploc = dfloc$pval 
    }else if (fdr_method=='BH'){
      ploc = p.adjust(dfloc$pval,method='BH')
    }else{
      stop('unknown FDR method')
    }
    num_TP = sum(ploc<=alpha & state == 'osc')
    num_FP = sum(ploc<=alpha & state == 'non_osc')
    return(cbind(methods[ind,],data.frame(TPR=num_TP/num_P,
                                 FPR=num_FP/num_N)))
  }) %>% rbindlist() %>% data.frame()
  return(ROCdf)
}