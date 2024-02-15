helper_auglattice_to_state=function(N1,N2,shift2,scale2){
  lat1=c(1:N1)/N1-1/N1
  lat2=c(1:N2)/N2-1/N1
  lat2=shift2+scale2*lat2
  mt=c(lat1,lat2)
  return(mt)
}