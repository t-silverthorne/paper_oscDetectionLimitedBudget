#' Helper function converts parameterization of two lattices into single lattice
convert_2lattice_to_state=function(shift1,shift2,scale1,scale2,lat1,lat2){
  x1 = shift1+scale1*lat1
  x2 = shift2+scale2*lat2
  x  = c(x1,x2) %% 1 # convention, keep measurement times inside study
  return(x)
}