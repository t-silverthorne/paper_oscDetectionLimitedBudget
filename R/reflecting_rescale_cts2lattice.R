reflecting_rescale_cts2lattice=function(shift,scale,rnd_scale,lat){
  if( (shift*(1+rnd_scale) > 0) & 
      (all(shift+scale*(1+rnd_scale)*lat>0))&
      (all(shift+scale*(1+rnd_scale)*lat<1))){
    scale_new = scale*(1+rnd_scale)
  }else if ( (shift*(1-rnd_scale) > 0) & 
             (all(shift+scale*(1-rnd_scale)*lat>0))&
             (all(shift+scale*(1-rnd_scale)*lat<1))){
    scale_new = scale*(1-rnd_scale)
  }else{
    scale_new=shift
  }
  return(scale_new)
}