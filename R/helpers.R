helper_auglattice_to_state=function(N1,N2,shift2,scale2){
  lat1=c(1:N1)/N1-1/N1
  lat2=c(1:N2)/N2-1/N1
  lat2=shift2+scale2*lat2
  mt=c(lat1,lat2)
  return(mt)
}

helper_unif_update_shift=function(scale,shift,lat,tscale){
  if (any(scale*lat + shift <0)  | any(scale*lat + shift >1)){
    stop('invalid initial state')
  }else{
    # TODO: replace with simpler parameterization 
    lmin=min(lat)
    lmax=max(lat) 
   
    sh_min = -scale*lmin
    sh_max = 1-scale*lmax
     
    eps_n = min(tscale/2,shift-sh_min)
    eps_p = min(tscale/2,sh_max-shift)
   
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    }
    shift_new = runif(1,shift-eps_n,shift+eps_p)
    return(shift_new)
  }
}

helper_unif_update_shift_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (scale*(Npts-1)+shift > Nfine | shift<0 | scale < 0){
    stop('invalid initial state')
  }else{
    sh_min = 0
    sh_max = Nfine-1 -(Npts-1)*scale
    
    xL = max(sh_min,shift-floor(tscale/2))
    xR = min(sh_max,shift+floor(tscale/2))
    
    shift_new = sample(c(xL:xR),1)
    return(shift_new)
  }
}


helper_unif_update_scale=function(scale,shift,lat,tscale){
  if (any(scale*lat + shift <0)  | any(scale*lat + shift >1)){
    stop('invalid initial state')
  }else{
    # TODO: replace with simpler parameterization 
    lmin=min(lat)
    lmax=max(lat) 
   
    sc_min = 0
    sc_max = (1-shift)/lmax
     
    eps_n = min(tscale/2,scale-sc_min)
    eps_p = min(tscale/2,sc_max-scale)
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    } 
    scale_new = runif(1,scale-eps_n,scale+eps_p)
    return(scale_new)
  }
}


helper_unif_update_scale_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (tscale<2){
    stop('for discrete problem, tscale must be >= 2')
  }
  if (scale*(Npts-1)+shift > Nfine | shift<0 ){
    stop('invalid initial state')
  }else{
    sc_min = 1 
    sc_max = floor((Nfine-shift-1)/(Npts-1)) #TODO: verify you don't need floor here
    
    xL = max(sc_min,scale-floor(tscale/2))
    xR = min(sc_max,scale+floor(tscale/2))
    
    scale_new = sample(c(xL:xR),1)
  }
  return(scale_new)
}

helper_gap_penalty=function(mt){
  mt=sort(unique(mt))
  if (any(mt>1) | any(mt<0)){
    stop('measurement times must be in [0,1]')
  }
  d1 = 1-(mt[length(mt)]-mt[1])
  return(max(c(diff(mt),d1)))
}

xinds_from_lat1lat2_pars=function(N1,dx1,N2,dx2,xshift2){
  lat1  = 1+dx1*c(0:(N1-1)) 
  lat2  = 1+xshift2 + dx2*c(0:(N2-1))
  xinds = c(lat1,lat2)
  return(xinds)
}

wrapper_sweep_lattice=function(optim_routine,Nvals,Nmeas,cts_flag,dvar0,...){
  fval_best = Inf
  res_best  = NULL
  N1_best   = NULL
  time_tot  = 0
  all_fvals = NULL 
  for (n1 in Nvals){
    n2=Nmeas-n1
    if (cts_flag){
      dvar0[['lat1']]= c(1:n1)/n1 - 1/n1
      dvar0[['lat2']]= c(1:n2)/n2 - 1/n2
      res_loc = optim_routine(dvar0,...)
    }else{
      dvar0$N1=n1
      dvar0$N2=n2
      res_loc = optim_routine(dvar0,...)
    }
    time_tot=time_tot+as.numeric(res_loc$runtime)
    all_fvals=rbind(all_fvals,data.frame(fval=res_loc$fvalue,N1=n1))
    if (res_loc$fvalue<fval_best){
      res_best  = res_loc
      N1_best   = n1
      fval_best = res_loc$fvalue
    }
  }
  return(list(res_best=res_best,N1_best=N1_best,
              time_tot=time_tot,
              all_fvals=all_fvals))
}

helper_extract_info_from_result=function(res,tag,Nmeas=gset$Nmeas,runtime_tot,cfunlabel='power',...){
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
  
  if (cfunlabel=='power'){
    return(annmatrix(
                 x = matrix(sort(res$mtvalue),nrow=1),
              rann = data.frame(power   = 1-res$fvalue,
                                runtime = runtime,
                                tag     = tag,
                                cts     = tagsp[[1]][1],
                                lattice = tagsp[[1]][2],
                                solver  = tagsp[[1]][3]),
              cann = data.frame(time=c(1:Nmeas))))
  }else if(cfunlabel=='ncp'){
    return(annmatrix(
                 x = matrix(sort(res$mtvalue),nrow=1),
              rann = data.frame(ncp     = -res$fvalue,
                                runtime = runtime,
                                tag     = tag,
                                cts     = tagsp[[1]][1],
                                lattice = tagsp[[1]][2],
                                solver  = tagsp[[1]][3]),
              cann = data.frame(time=c(1:Nmeas))))
  }
}

helper_update_amm=function(amm,Lnow,gset,tag,wrapped=F,...){
  for (ii in c(1:gset$nrep)){
    res = Lnow[[ii]]
    if (wrapped){
      amm = rbind(amm,helper_extract_info_from_result(res=res$res_best,tag=tag,
                                               runtime_tot=res$time_tot,...)) 
    }else{
      amm = rbind(amm,helper_extract_info_from_result(res=res,tag=tag,...)) 
    }
  }
  return(amm)  
}


