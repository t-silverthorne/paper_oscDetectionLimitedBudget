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