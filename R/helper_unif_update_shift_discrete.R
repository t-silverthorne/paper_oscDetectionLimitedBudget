helper_unif_update_shift_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (scale*(Npts-1)+shift > Nfine | shift<0 | scale < 0){
    stop('invalid initial state')
  }else{
     
    eps_n = min(floor(tscale/2),shift)
    eps_p = min(floor(tscale/2),floor((Nfine-1-shift)/(Npts-1)/scale))
   
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    }
    shift_new = sample(c((shift-eps_n):(shift+eps_p)),1)
    return(shift_new)
  }
}