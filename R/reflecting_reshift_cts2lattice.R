reflecting_reshift_cts2lattice=function(shift,scale,rnd_shift,lat){
  shift_prop  = shift + rnd_shift 
  shift_propn = shift - rnd_shift
  lat_prop    = shift_prop  + scale*lat 
  lat_propn   = shift_propn + scale*lat
  if (all( (lat_prop<=1) & (lat_prop>=0))){
    shift=shift_prop
  }else if(all( (lat_propn<=1) & (lat_propn>=0))){
    shift=shift_propn
  }else{
    stop('unable to transition')
  }
  return(shift) 
}
  