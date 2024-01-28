reflecting_reshift_cts2lattice=function(shift,scale,rnd_shift,lat){
  if (all(shift+rnd_shift+scale*lat >0) & all(shift+rnd_shift+scale*lat <1)){
    shift_new=shift+rnd_shift
  }else if(all(shift-rnd_shift+scale*lat >0) & all(shift-rnd_shift+scale*lat <1)){
    shift_new=shift-rnd_shift
  }else{
    shift_new=shift
  }
  return(shift_new)
}