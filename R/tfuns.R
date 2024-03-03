tfun_auglattice=function(N1,N2,shift2,scale2,freqs,control){
  if (N1<control$N1min){
    stop('too few points in lattice1')
  }
  if (N1>control$N1max){
    stop('too many points in lattice1')
  }
  # update number of points in each lattice
  N1cands=c(N1,min(N1+5,control$N1max),max(N1-5,control$N1min))
  N1new = sample(N1cands,1) 
  N2new = N2+N1-N1new
  
  shift2new = runif(1,0,1/N1new)
  scale2new = runif(1,0,N2new*(1-shift2new)/(N2new-1))
  return(c(N1new,N2new,shift2new,scale2new))
}

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
}else{
  stop('unknown control$tfun_choice')
  
}
}

#' requires control$tfun_choice
#'          control$tfun_mean
#'          control$tfun_sd 
tfun_svdpower=function(mt,control){
  if(control$tfun_choice=='brownian-torus'){
    mt=mt+rnorm(length(mt),mean=control$tfun_mean,sd = control$tfun_sd) 
    mt=mt%%1
  }else{
    stop('unknown control$tfun_choice')
  }
  return(mt)
}


#' requires control$tfun_choice
tfun_svdpower_discrete=function(xinds,Nfine,control){
  if (control$tfun_choice=='single-flip'){
    flip_index        = sample(length(xinds),1)
    new_val           = sample(setdiff(c(1:Nfine),xinds),1) 
    xinds[flip_index] = new_val
    #TODO: add discretized Brownian motion
  }else{
    stop('unknown control$tfun_choice')
  }
  return(xinds)
}

#' requires control$tfun_choice
tfun_2lattice_svdpower_discrete=function(N1,dx1,N2,dx2,xshift2,Nfine,control){
  if(control$tfun_choice=='unif-with-bdry-discrete'){
    accept_state=F
    while(accept_state==F){
      dx1_new     =  helper_unif_update_scale_discrete(dx1,0,N1,control$tscale,Nfine)
      dx2_new     =  helper_unif_update_scale_discrete(dx2,xshift2,N2,control$tscale,Nfine)
      xshift2_new =  helper_unif_update_shift_discrete(dx2_new,xshift2,N2,control$tscale,Nfine)
      
      if(!anyDuplicated(xinds_from_lat1lat2_pars(N1,dx1_new,N2,dx2_new,xshift2_new))){
        accept_state = T
      }
    }
  }else if (control$tfun_choice=='reflecting-brownian'){
    stop('reflecting-brownian is depracated')
  }else{
    stop('unknown control$tfun_choice')
  }
  return(list(dx1=dx1_new,dx2=dx2_new,xshift2=xshift2_new)) 
  
}