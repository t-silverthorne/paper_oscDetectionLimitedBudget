#' Requires control$tfun_choice
#'          control$mean_scale1
#'          control$sd_scale1
#'          control$mean_scale2
#'          control$sd_scale2
#'          control$mean_shift2
#'          control$sd_shift2
tfun_2lattice_svdpower=function(shift1,shift2,scale1,scale2,lat1,lat2,
                                control){
if(control$tfun_choice=='unif-with-bdry'){
  scale1=helper_unif_update_scale(scale1,shift1,lat1,control$tscale_unif_with_bdry)
  scale2=helper_unif_update_scale(scale2,shift2,lat2,control$tscale_unif_with_bdry)
  shift2=helper_unif_update_shift(scale2,shift2,lat2,control$tscale_unif_with_bdry)
  return(list(scale1=scale1,
              scale2=scale2,
              shift2=shift2)) 
}else if (control$tfun_choice=='default'){ # keeps shift1=0
  stop('Reflecting boundaries not yet tested/implemented')
  # generate random numbers
  rnd_scale1 = rnorm(1,control$mean_scale1,control$sd_scale1)  
  rnd_scale2 = rnorm(1,control$mean_scale2,control$sd_scale2)
  rnd_shift2 = rnorm(1,control$mean_shift2,control$sd_shift2)
 
  scale1     = reflecting_rescale_cts2lattice(shift1,scale1,rnd_scale1,lat1)
  scale2     = reflecting_rescale_cts2lattice(shift2,scale2,rnd_scale2,lat2) 
  shift2     = reflecting_reshift_cts2lattice(shift2,scale2,rnd_shift2,lat2)
  
  return(list(scale1=scale1,
              scale2=scale2,
              shift2=shift2)) 
}else{
  stop('unknown control$tfun_choice')
  
}
}