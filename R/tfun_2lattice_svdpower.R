#' Requires control$tfun_choice
#'          control$mean_scale1
#'          control$sd_scale1
#'          control$mean_scale2
#'          control$sd_scale2
#'          control$mean_shift2
#'          control$sd_shift2
tfun_2lattice_svdpower=function(shift1,shift2,scale1,scale2,lat1,lat2,
                                control){
if (control$tfun_choice=='default'){ # keeps shift1=0
  rnd_scale1 = rnorm(1,control$mean_scale1,control$sd_scale1)  
  rnd_scale2 = rnorm(1,control$mean_scale2,control$sd_scale2)
  rnd_shift2 = rnorm(1,control$mean_shift2,control$sd_shift2)
 
  scale1_new = reflecting_rescale_cts2lattice(shift1,scale1,rnd_scale1,lat1) 
  scale2_new = reflecting_rescale_cts2lattice(shift2,scale2,rnd_scale2,lat2) 
  shift2_new = reflecting_reshift_cts2lattice(shift2,scale2_new,rnd_shift2,lat2) #important, output of prev influences this test 
  
  return(list(scale1=scale1_new,
              scale2=scale2_new,
              shift2=shift2_new)) 
}else{
  stop('unknown control$tfun_choice')
  
}
}