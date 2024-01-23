#' Helper function for continuous power optimization
#' @export
costfun_svdpower=function(mt,freqs,Amp=1,alpha=.05){
  N   = length(mt)
  f0  = qf(p=1-alpha,df1=2,df2=N-3)
  
  return(freqs %>% sapply(function(freq){
  getReducedFIM(mt,list(freq=freq)) %>% {Amp^2*getMinEig(.)} %>% 
                  {1 - pf(q=f0,df1=2,df2=N-3,ncp=.)} 
  })%>% min())
}