runCVXtest <- function(n){
  t=runif(n) %>% matrix(nrow=n)
  s=runif(n) %>% matrix(nrow=n)
  q=runif(1)
  
  cfun = function(t,p){
    cvec = cos(2*pi*p$freq*t)
    svec = sin(2*pi*p$freq*t)
    
    ( t(cvec)%*%cvec )*( t(svec)%*%svec) - (t(cvec)%*%svec)^2
  }
  
  param=list(freq=runif(1)*10)
  return(q*cfun(t,param) + (1-q)*cfun(s,param) - cfun(q*t + (1-q)*s,param))
}

n=10
N=1e3
replicate(N,{runCVXtest(n)}) %>% summary()