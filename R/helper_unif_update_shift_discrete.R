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