getPcQuadCsts=function(fmin,fmax,Nfreq,
                       w_reg,w_wc,Nfine){
  QC=list() # list of quadratic constraints
 
  # matrices for regularization term 
  tau  = c(1:Nfine)/Nfine -1/Nfine 
  dtau_mat          = outer(tau,tau,'-')
  tau_prod_mat      = outer(tau,tau,'*')
  ncp_reg_mat       = ((sin(4*pi*fmax*dtau_mat) - sin(4*pi*fmin*dtau_mat))*(4*tau_prod_mat*pi^2 + 1))/(4*dtau_mat)
  diag(ncp_reg_mat) = pi*(4*pi^2*tau^2 + 1)*(fmax - fmin) 
  
  if (w_wc==0){
    qcmat  = w_reg*ncp_reg_mat
    qcmat  = cbind(rbind(qcmat,rep(0,Nfine)),rep(0,Nfine+1))
    qc     = c(rep(0,Nfine),-1)
    beta   = 0
    sense  = '<'
    QC[[1]] = list(Qc=qcmat,q=qc,rhs=beta,sense=sense)
  }else if (w_wc>0){
    fvec = seq(from=fmin,to=fmax,length.out=Nfreq)
    for (ii in c(1:Nfreq)){
      freq  = fvec[ii]
      cvec  = matrix(cos(2*pi*freq*tau),nrow=Nfine)
      svec  = matrix(sin(2*pi*freq*tau),nrow=Nfine)
      
      a11   = cvec*cvec
      a22   = svec*svec
      a12   = cvec*svec
      
      qcmat  = w_wc*(a11%*%t(a11) + a22%*%t(a22) +4*a12%*%t(a12)-a11%*%t(a22) - a22%*%t(a11))+w_reg*ncp_reg_mat
      qcmat  = cbind(rbind(qcmat,rep(0,Nfine)),rep(0,Nfine+1))
      qc     = c(rep(0,Nfine),-1)
      beta   = 0
      sense  = '<'
      QC[[ii]] = list(Qc=qcmat,q=qc,rhs=beta,sense=sense)
    }
    return(QC) 
  }
}