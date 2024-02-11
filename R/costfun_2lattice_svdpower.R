#' Helper function for continuous power optimization, restricted to two lattice subspace
#' @export
costfun_2lattice_svdpower=function(shift1,shift2,scale1,scale2,lat1,lat2,
                                   freqs,Amp=1,alpha=.05,regL1=0,...){
  mt=convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2)
  return(costfun_svdpower(mt,freqs,Amp,alpha,regL1,...))
}