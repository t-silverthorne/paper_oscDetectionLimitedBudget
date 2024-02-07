helper_unif_update_scale_discrete=function(scale,shift,Npts,tscale,Nfine){
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