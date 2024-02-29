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