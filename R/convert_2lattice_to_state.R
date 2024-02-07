#' Helper function converts parameterization of two lattices into single lattice
convert_2lattice_to_state=function(shift1,shift2,scale1,scale2,lat1,lat2){
  x1 = shift1+scale1*lat1
  x2 = shift2+scale2*lat2
  # keep measurement times inside study, necessary for L-BFGS-B solver but not 
  # simulated annealing (in simulated annealing, times stay inside study because
  # of design of transition function) 
  x  = c(x1,x2) 
  return(x)
}