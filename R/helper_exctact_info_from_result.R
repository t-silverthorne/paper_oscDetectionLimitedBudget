helper_extract_info_from_result=function(res,tag,Nmeas=gset$Nmeas,runtime_tot,...){
  tagsp=strsplit(tag,'_')
  if(hasArg('runtime_tot')){
    # lubridate ensures that units are always in seconds
    runtime=runtime_tot %>% lubridate::as.duration() %>% as.numeric()
  }else{
    runtime=res$runtime %>% lubridate::as.duration() %>% as.numeric()
  }
  if (length(res$mtvalue)<Nmeas){
    stop('degenerate measurement times')
  }
  return(annmatrix(
               x = matrix(sort(res$mtvalue),nrow=1),
            rann = data.frame(power   = 1-res$fvalue,
                              runtime = runtime,
                              tag     = tag,
                              cts     = tagsp[[1]][1],
                              lattice = tagsp[[1]][2],
                              solver  = tagsp[[1]][3]),
            cann = data.frame(time=c(1:Nmeas))))
}