xinds_from_lat1lat2_pars=function(N1,dx1,N2,dx2){
  lat1  = dx1*c(0:(N1-1)) 
  lat2  = xshift2 + dx2*c(0:(N2-1))
  xinds = c(lat1,lat2)
}