getPower <- function(t,param,alpha=.05){
# return power of one-frequency cosinor model
  Amp    = param[['Amp']];
  freq   = param[['freq']];
  acro   = param[['acro']];
  N      = length(t)
  
  cvec   = Amp*cos(2*pi*freq*t-acro);
  lambda = as.numeric(t(cvec)%*%cvec)
  
  f0     = qf(p=1-alpha,df1=2,df2=N-3)
  return(1 - pf(q=f0,df1=2,df2=N-3,ncp=lambda))
}
getPower_acrovec <- function(t,param,acrovec,alpha=.05){
# loop previous
  powervec = rep(NaN,length(acrovec))
  for (ii in 1:length(acrovec)){
    param$acro   = acrovec[ii]
    powervec[ii] = getPower(t,param,alpha)
  }
  return(powervec)
}

getMinPower <- function(t,param,alpha=.05){
# find minimum power over the interval from [0,2pi] 
  return(pracma::fminbnd(
        f = function(phi){param['acro']=phi
            return(getPower(t,param,alpha))},
        a = 0,
        b = 2*pi) %>% {.$fmin})
}

optPowerDirect<-function(t,param,alpha=.05){
  return(pracma::fmincon(x0=t,
                  fn=function(t){-getMinPower(t,param,alpha)},
                  lb=rep(0,length(t)),
                  ub=rep(1,length(t))))
}

optPowerLambda <- function(t,param,alpha=.05){
  N = length(t)
  pracma::fmincon(x0=t,
                  fn=function(t){
                      t     = matrix(t,nrow=N)
                      cvec  = cos(2*pi*param$freq*t)
                      svec  = sin(2*pi*param$freq*t)
                      fdet  = (t(cvec)%*%cvec)*(t(svec)%*%svec) - (t(cvec)%*%svec)^2
                      return(-0.5*(N - sqrt(N^2 - 4*fdet)))
                  },
                  lb=rep(0,length(t)),
                  ub=rep(1,length(t)))
}

getFIM <- function(t,param){
  Amp    = param[['Amp']];
  freq   = param[['freq']];
  acro   = param[['acro']];
  N      = length(t)
  
  Y = matrix(c(rep(1,length(t)),
           cos(2*pi*freq*t),
           sin(2*pi*freq*t)),
           nrow=3,byrow=T)

  Y %*% t(Y)
}
getReducedFIM <- function(t,param){
  Amp    = param[['Amp']];
  freq   = param[['freq']];
  acro   = param[['acro']];
  N      = length(t)
  
  Y = matrix(c(cos(2*pi*freq*t),
           sin(2*pi*freq*t)),
           nrow=2,byrow=T)
  Y %*% t(Y)
}
