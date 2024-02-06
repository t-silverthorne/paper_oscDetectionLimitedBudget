helper_unif_update_scale_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (scale*(Npts-1)+shift > Nfine | shift<0 ){
    stop('invalid initial state')
  }else{
    sc_min = 1 
    sc_max = (Nfine-shift-1)/(Npts-1) #TODO: verify you don't need floor here
     
    eps_n = min(floor(tscale/2),floor(scale-sc_min))
    eps_p = min(floor(tscale/2),floor(sc_max-scale))
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    } 
    scale_new = sample(c((scale-eps_n):(scale+eps_p)),1)
  }
  return(scale_new)
}