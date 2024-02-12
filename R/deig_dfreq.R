deig_dfreq=function(mt,freq){
  cvec=cos(2*pi*freq*mt)
  svec=sin(2*pi*freq*mt)
  
  M11=t(cvec)%*%cvec
  M22=t(svec)%*%svec
  M12=t(cvec)%*%svec
  
  M  = matrix(c(M11,M12,M12,M22),nrow=2)
  M
  eigs = eigen(M,symmetric = T)
  vmin = eigs$vectors[,which.min(eigs$values)]
  
  dM22     = sum(4*pi*mt*cvec*svec)
  dM12     = sum(2*pi*mt*(2*cvec*cvec-1))
  dMdfreq = matrix(c(-dM22,dM12,dM12,dM22),nrow=2,byrow=T)
  
  return(t(vmin)%*%dMdfreq%*%vmin/(t(vmin)%*%vmin))
}



