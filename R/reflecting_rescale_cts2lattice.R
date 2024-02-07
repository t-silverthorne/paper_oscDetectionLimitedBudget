reflecting_rescale_cts2lattice=function(shift,scale,rnd_scale,lat){
  scale_prop  = scale+rnd_scale
  scale_propn = scale-rnd_scale
  lat_prop    = shift + scale_prop*lat 
  lat_propn   = shift + scale_propn*lat
  if (all( (lat_prop<=1) & (lat_prop>=0))){
    scale=scale_prop
  }else if(all( (lat_propn<=1) & (lat_propn>=0))){
    stop('negative')
    scale=scale_propn
  }else{
    stop('unable to transition')
  }
  return(scale) 
}
  