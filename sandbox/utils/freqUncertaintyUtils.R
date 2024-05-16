minLambdaOverFreq = function(tvec,param,fmin,fmax,method,delf=.01){
  minLambda=NaN
  if (strcmp(method,'fast')){
      minLambda = pracma::fminbnd(
        function(f){
          param$freq = f
          return(getMinEig(tvec,param))
          },
      a=fmin,b=fmax) %>% {.$fmin}
  } else if (strcmp(method,'slow')){
      freqlistloc=seq(from=fmin,to=fmax,by=delf)
      minLambda = freqlistloc %>% sapply(function(f){
        param$freq=f 
        return(getMinEig(tvec,param))}
      ) %>% min()
  }
  return(minLambda)
}
minPowerOverFreq = function(tvec,param,fmin,fmax,method,delf=.01){
  minLambda=NaN
  if (strcmp(method,'fast')){
      minLambda = pracma::fminbnd(
        function(f){param$freq = f
        return(getMinPower(tvec,param))},
      a=fmin,b=fmax) %>% {.$fmin}
  } else if (strcmp(method,'slow')){
      freqlistloc=seq(from=fmin,to=fmax,by=delf)
      minLambda = freqlist %>% sapply(function(f){
        param$freq=f 
        return(getMinPower(tvec,param))}
      ) %>% min()
  }
  return(minLambda)
}