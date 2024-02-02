helper_unif_update_shift_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (scale*(Npts-1)+shift > Nfine | shift<0 ){
    stop('invalid initial state')
  }else{
    lmin=shift
    lmax=(Npts-1)
   
    sh_min = -scale*lmin
    sh_max = Nfine-scale*lmax # TODO: verify you don't need floor here
     
    eps_n = floor(min(tscale/2,shift-sh_min))
    eps_p = floor(min(tscale/2,sh_max-shift))
   
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    }
    shift_new = sample(c((shift-eps_n):(shift+eps_p)),1)
    return(shift_new)
  }
}