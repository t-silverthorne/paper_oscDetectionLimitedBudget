#' requires control$tfun_choice
tfun_2lattice_svdpower_discrete=function(N1,dx1,N2,dx2,xshift2,Nfine,control){
  if (control$tfun_choice=='reflecting-brownian'){
    accept_state=F
    while (accept_state==F){
      ddx1    = sample(c(-1,0,1),1)
      ddx2    = sample(c(-1,0,1),1)
      dshift2 = sample(c(-1,0,1),1) 
    
      if (dx1+ddx1 < 1 |(dx1+ddx1)*(N1-1)>Nfine){
        dx1_new = dx1+ddx1 
      }else{
        dx1_new = dx1-ddx1 
      }
      
      if (dx2+ddx2 < 1 |xshift2 + (dx2+ddx2)*(N2-1)>Nfine){
        dx2_new = dx2+ddx2 
      }else{
        dx2_new = dx2-ddx2  #TODO: check this is exhaustive
      }
    
      if( xshift2 + dshift2>0 &(xshift2+dshift2 + (dx2_new)*(N2-1)<=Nfine)){
        xshift2_new = xshift2 + dshift2
      }else if( xshift2 - dshift2>0 &(xshift2-dshift2 + (dx2_new)*(N2-1)<=Nfine)){
        xshift2_new = xshift2 - dshift2
      }else{
        xshift2_new = xshift2  
      } 
      if(!anyDuplicated(xinds_from_lat1lat2_pars(N1,dx1_new,N2,dx2_new,xshift2_new))){
        accept_state = T
      }
    }
    return(list(dx1=dx1_new,dx2=dx2_new,xshift2=xshift2_new)) 
  }else{
    stop('unknown control$tfun_choice')
  }
  
}