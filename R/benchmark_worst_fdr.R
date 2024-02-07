benchmark_worst_fdr=function(mt,Amp,f0,Nmc,Nacro,Nfreq,fmin,fmax){
    f0    = 1 
    fmin  = 1
    Amp   = 2
    fmax  = 23.5
    Nfreq = 2^7
    Nmc   = 1e2
    freqs = seq(from=fmin,to=fmax,length.out=Nfreq)
    acros = seq(from=0,to=2*pi,length.out=2^3)
    
    acros %>%sapply(function(acro){
      X=t(replicate(Nmc,{Amp*cos(2*pi*f0*mt-acro)}))+
        matrix(rnorm(Nmc*length(mt)),nrow=Nmc)
     
      Lfit = freqs %>%as.list() %>% lapply(function(freq){
        rowCosinor(X,zts=mt,per=1/freq) %>% {data.frame(pvalue=t(.$pvalue))}})
     
      pmat=Lfit%>% rbindlist() %>% as.matrix() %>% t()
      qmat=matrix(replicate(Nmc*Nfreq,{NaN}),nrow=Nmc)
      for (ii in c(1:Nmc)){
        qmat[ii,]=p.adjust(pmat[ii,],'BY')
      }
      is_osc = apply(qmat,1,function(x){any(x<.05)}) 
      if (length(is_osc)==Nmc){
        return(mean(is_osc))
      }else{
        stop('inconsistent dimension, did you transpose something?')
      }
    }) %>% min()
}