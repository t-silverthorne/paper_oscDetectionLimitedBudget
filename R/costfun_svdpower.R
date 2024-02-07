#' Helper function for continuous power optimization
#' @export
costfun_svdpower=function(mt,freqs,Amp=1,alpha=.05,regL1=0){
  N   = length(mt)
  f0  = qf(p=1-alpha,df1=2,df2=N-3)
  
  pwrs = freqs %>% sapply(function(freq){
  getReducedFIM(mt,list(freq=freq)) %>% {Amp^2*getMinEig(.,is_symmetric=T)} %>% 
                  {1 - pf(q=f0,df1=2,df2=N-3,ncp=.)} 
  })
  if (regL1>0){
    cval = min(pwrs)+regL1*mean(pwrs)  
  }else{
    cval = min(pwrs)
  }
  return(cval)
}