
rowCosinor <- function(theData, zts, per=24) {
  Y <- as.matrix(theData)
  
  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)
  
  betas <- solve(t(X) %*% X,t(X) %*% t(Y))
  
  phases     <- atan2(betas[2,], betas[3,]) %% (2*pi)
  amplitudes <- sqrt(betas[2,]*betas[2,] + betas[3,]*betas[3,])
  
  fits <- t(X %*% betas)
  
  SStot <- rowSums((Y - rowMeans(Y))^2)
  SSres <- rowSums((fits-Y)^2)
  Rsqs  <- 1 - (SSres/SStot)
  
  SSmod <- SStot - SSres
  DFres <- ncol(theData) - 3
  DFmod <- 2
  MSres <- SSres / DFres
  MSmod <- SSmod / DFmod
  Fstatistic <- MSmod / MSres
  
  pval <- pf(Fstatistic, DFmod, DFres, lower.tail=FALSE)
  
  data.frame(phase=phases, amplitude=amplitudes, mesor=betas[1,],
             rsq=Rsqs, statistic=Fstatistic, pvalue=pval
  )
}

sols = readRDS('results/data/MCperiodogram/hiresSols.RDS')
tvec=sols[sols@wreg==1 & sols@drts==Inf & sols@Nmeas==32,]
tvec=tvec[tvec<Inf]
#tvec=c(1:Nmeas)/Nmeas-1/Nmeas
freq=3.01
Amp=.75

evalWorstPower(tvec,freq,Amp)

Nmc=1e4
p_osc=.5
Nacro=2^5
Nmeas=length(tvec)
acro=seq(0,2*pi,length.out=Nacro+1)
acro=acro[1:Nacro]
pwrvec=sapply(c(1:Nacro),function(ii){
  param=list(freq=freq,acro=acro[ii],Amp=Amp)
  evalExactPower(tvec,param) 
})

TPRs=sapply(c(1:Nacro),function(ii){
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  Ydat[state=='osc',] =Ydat[state=='osc',] +
    t(replicate(sum(state=='osc'),Amp*cos(2*pi*tvec*freq-acro[ii])))

   pvdf = data.frame(pval=rowCosinor(Ydat,tvec,per=1/freq) %>% {.$pvalue},
                     state=state)
   pval   = pvdf$pval
   ostate = pvdf$state
   
   # record FPR and TPR at alpha=.05
   num_P   = sum(ostate=='osc') 
   num_N   = sum(ostate=='non_osc') 
   num_TP  = sum(ostate=='osc'     & pval < .05) 
   num_FP  = sum(ostate=='non_osc' & pval < .05)
   TPR     = num_TP/num_P
   return(TPR)
})
mean(TPRs)
mean(pwrvec)

min(pwrvec)
min(TPRs)

df1=data.frame(acro=acro,fom=pwrvec,type='power')
df2=data.frame(acro=acro,fom=TPRs,type='TPR')
df1
df2
rbind(df1,df2) %>% ggplot(aes(x=acro,y=fom,group=type,color=type))+geom_line()+ylim(c(0,1))
