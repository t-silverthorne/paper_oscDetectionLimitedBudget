#' Get p-values from Lomb Scargle periodogram 
#' @param alpha_vals values of threshold parameter
#' @param Nmc number of signals to include in dataset
#' @param tvec vector of measurement times
#' @param Amp amplitude of signal
#' @param freq frequency of signal
#' @param p_osc what proportion of the simulated data contains oscillations
#' @param mc_cores number of cores to be used for parallelization
#' @return state vector indicating which genes are osc and which are not osc
#' a dataframe of simulated p-values 
#' 
#' @author Turner Silverthorne
#' @export
simLombScarglePvals=function(alpha_vals,Nmc,tvec,
                       Amp=1,freq=1,p_osc=.1,mc_cores=1){
  # generate simulated data
  Nmeas               = length(tvec)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*tvec,'-'))
  
  # compute p-values from Lomb-Scargle
  pvdf = c(1:dim(Ydat)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=tvec,plot=F,normalize = 'standard')
    lomb_press     = lsp(x,times=tvec,plot=F,normalize = 'press')
    return(rbind(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]),
                 data.frame(p_method='press',pval =lomb_press$p.value,state=state[ii])))
  }) %>% rbindlist() %>% data.frame()
  return(list(state=state,pvdf=pvdf))
}