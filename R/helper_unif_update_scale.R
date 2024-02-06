helper_unif_update_scale=function(scale,shift,lat,tscale){
  if (any(scale*lat + shift <0)  | any(scale*lat + shift >1)){
    stop('invalid initial state')
  }else{
    #TODO: check parameterization is consistent
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