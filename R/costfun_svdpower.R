#' Helper function for continuous power optimization
#' @export
costfun_svdpower=function(mt,freqs,Amp=1,alpha=.05,regL1=0,cfuntype='power'){
  N   = length(mt)
  f0  = qf(p=1-alpha,df1=2,df2=N-3)

  # depending on cfuntype, return non-centrality parameter or power
  if (cfuntype=='ncp'){
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {getMinEig(.,is_symmetric=T)} 
    })
    
  } else if (cfuntype=='power'){
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {Amp^2*getMinEig(.,is_symmetric=T)} %>% 
                    {1 - pf(q=f0,df1=2,df2=N-3,ncp=.)} 
    })
  }else{
    stop('unknown cfuntype')
  }
  if (regL1>0){
    cval = min(cvals)+regL1*mean(cvals)  
  }else{
    cval = min(cvals)
  }
  return(cval)
}