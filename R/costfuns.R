#'Cost function for continuous-time power optimisation
#' @param mt vector of measurement times 
#' @param freqs vector of frequencies at which cost function should be evaluated
#' @param Amp amplitude of signal, only relevant if \code{cfuntype=='power'}
#' @param alpha type I error cutoff, only relevant if \code{cfuntype=='power'}
#' @param weight_ncp factor for scaling non-centrality parameter in multi-objective scalarization. 
#' @param regL1 factor for scaling L1-frequency regularisation in multi-objective scalarization. 
#' @param regFder factor for scaling frequency derivative in multi-objective scalarization. 
#' @param gapPenalty factor for scaling adding a penalty for gaps between consecutive measurement times
#' @param leveragePenalty currently depracated, a factor for penalising by statistical leverage
#' @param cfuntype either 'power' or 'ncp' if non-centrality parameter is to be optimised directly
#' 
#' @return minimum of cost function evaluated at \code{freqs}
#' @author Turner Silverthorne
#' @export
costfun_svdpower=function(mt,freqs,
                          Amp=1,alpha=.05,weight_ncp=1,regL1=0,regFder=0,
                          gapPenalty=0,leveragePenalty=0,
                          cfuntype='power'){
  N   = length(mt)
  if (regFder!=0 & cfuntype=='power'){
    stop('Use cfuntype=ncp, frequency regularization for power has not been implemented.')
  }
  if (cfuntype=='ncp'){
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {getMinEig(.,is_symmetric=T)} 
    })
    
  } else if (cfuntype=='power'){
    f0  = qf(p=1-alpha,df1=2,df2=N-3)
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {Amp^2*getMinEig(.,is_symmetric=T)} %>% 
                    {1 - pf(q=f0,df1=2,df2=N-3,ncp=.)} 
    })
  }else{
    stop('unknown cfuntype')
  }
  cval = weight_ncp*min(cvals)
  if (regL1>0){
    cval = cval+regL1*mean(cvals)  
  }
  
  if (regFder!=0){
    dlambda_dfreq = freqs %>% sapply(function(freq){getDeigDFreq(mt,freq)}) %>% {.^2} %>% sum()
    dlambda_dfreq = sqrt(dlambda_dfreq/length(freqs))
    if (is.numeric(regFder)&regFder>0){
      cval = cval-regFder*dlambda_dfreq
    }
    else if(regFder=='2norm'){
      cval = dlambda_dfreq
    }else{
      stop('unrecognised regFder setting')
    }
  }
  
  if (gapPenalty>0){
    cval=cval-gapPenalty*helper_gap_penalty(mt) 
  }
  if (leveragePenalty>0){
    stop('leverage depracated')
  }
  return(cval)
}

