#' Helper function for continuous power optimization
#' @export
costfun_svdpower=function(mt,freqs,Amp=1,alpha=.05,weight_ncp=1,regL1=0,regFder=0,gapPenalty=0,leveragePenalty=0,
                          cfuntype='power'){
  N   = length(mt)
  f0  = qf(p=1-alpha,df1=2,df2=N-3)
  if (regFder!=0 & cfuntype=='power'){
    stop('Use cfuntype=ncp, frequency regularization for power has not been implemented.')
  }
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
  cval = weight_ncp*min(cvals)
  if (regL1>0){
    cval = cval+regL1*mean(cvals)  
  }
  
  if (regFder>0){
    #dlambda_dfreq = freqs %>% sapply(function(freq){deig_dfreq(mt,freq)}) %>% abs() %>% max() 
    dlambda_dfreq = freqs %>% sapply(function(freq){deig_dfreq(mt,freq)}) %>% {.^2} %>% sum()
    dlambda_dfreq = sqrt(dlambda_dfreq/length(freqs))
    cval = cval-regFder*dlambda_dfreq
  }
  
  if (gapPenalty>0){
    cval=cval-gapPenalty*helper_gap_penalty(mt) 
  }
  if (leveragePenalty>0){
    maxlev = freqs %>% sapply(function(freq){
       cosinor_design_matrix(mt,freq) %>% {stats::hat(.,intercept=F)} %>% abs() %>% max()
    }) %>% max()
    cval =cval - leveragePenalty*maxlev
  }
  return(cval)
}