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

getMinEig <- function(t,param){
  R = getReducedFIM(t,param)
  return(eigen(R) %>% {.$values} %>% min())
}
multipwr=function(t,freqlist,param){
  freqlist %>% sapply(function(f){param$freq=f
  return(getMinPower(t,param))}) %>% min()
}

getMinEig_overFreq = function(t,freqlist,param){
  freqlist %>% sapply(function(f){param$freq=f 
  return(getMinEig(t,param))}
  ) %>% min()
}

getMeanEig_overFreq = function(t,freqlist,param){
  freqlist %>% sapply(function(f){param$freq=f 
  return(getMinEig(t,param))}
  ) %>% mean()
}


diffLambdaMin = function(t,param){
  # get eigenvector corresp. to minimal eigenvalue 
  M    = getReducedFIM(t,param)
  eM   = eigen(M)
  mind = which.min(eM$values)
  yv   = eM$vectors[,mind] 
  
  # compute entries of derivative mat 
  freq = param$freq
  om   = 2*pi*freq
  dA11 = -2*t(cos(om*t)) %*% (sin(om*t)*t)
  dA22 = 2*t(cos(om*t)) %*% (sin(om*t)*t)
  dA12 = t(cos(om*t)) %*% (cos(om*t)*t) - t(sin(om*t)) %*% (sin(om*t)*t)
  
  # use eigenvalue perturbation result
  dA   = 2*pi*matrix(c(dA11,dA12,dA12,dA22),nrow=2)
  
  return(t(yv)%*%dA%*%yv/(t(yv)%*%yv))
}

diffLambdaMin_overFreq = function(t,freqlist,param){
  der = freqlist %>% sapply(function(f){param$freq=f 
    return(diffLambdaMin(t,param))}
    )
  return(sqrt(t(der)%*%der))
}

getMinEig_overFreq_Reg = function(t,freqlist,param){
  minEig = freqlist %>% sapply(function(f){param$freq=f 
    return(getMinEig(t,param))}
    ) %>% min()
  reg    = diffLambdaMin_overFreq(t,freqlist,param)
  return(minEig-reg*param$regStrength)
}

getMinEigVariation = function(t,freqlist,param){
  return(-1*diffLambdaMin_overFreq(t,freqlist,param))
}

construct_poly_design       = function(delta,tvec1,tvec2){
  c((tvec1+c(delta))%%1,tvec2) %>% sort()
}
construct_sequential_design = function(tau,tvec1,tvec2){
  svec1 = tvec1*c(tau)
  svec2 = tvec2*(1-c(tau))
  return(c(svec1,c(tau)+svec2) %>% sort())
}

Jfun_delta = function(delta,tvec1,tvec2,freqlist,param){
  # cost fun for constructing polyrhythmic designs
  tvec  = construct_poly_design(delta,tvec1,tvec2) 
  if(param$method=='min'){
   J=getMinEig_overFreq(tvec,freqlist,param) 
  }else if(param$method=='min-reg'){
   J=getMinEig_overFreq_Reg(tvec,freqlist,param)
  }else if(param$method=='reg-only'){
   J=getMinEigVariation(tvec,freqlist,param)
  }else if(param$method=='mean'){
   J=getMeanEig_overFreq(tvec,freqlist,param)
  }
  return(J)
}

Jfun_tau = function(tau,tvec1,tvec2,freqlist,param){
  # cost fun for constructing sequential designs
  tvec = construct_sequential_design(tau,tvec1,tvec2)
  if(param$method=='min'){
   J=getMinEig_overFreq(tvec,freqlist,param) 
  }else if(param$method=='min-reg'){
   J=getMinEig_overFreq_Reg(tvec,freqlist,param)
  }else if(param$method=='reg-only'){
   J=getMinEigVariation(tvec,freqlist,param)
  }else if(param$method=='mean'){
   J=getMeanEig_overFreq(tvec,freqlist,param)
  }
  return(J)
}