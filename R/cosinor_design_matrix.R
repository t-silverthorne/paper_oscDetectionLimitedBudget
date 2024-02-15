cosinor_design_matrix=function(mt,freq){
  mt=as.matrix(mt,nrow=length(mt))
  return(cbind(matrix(rep(1,10),nrow=length(mt)),cos(2*pi*freq*mt),sin(2*pi*freq*mt)))
}