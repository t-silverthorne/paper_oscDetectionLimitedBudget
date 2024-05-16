evalAcroFreqGradL2 = function(tvec,prior='uniform',pars=list(fmin=1,fmax=24),Amp=1){
  if (prior=='uniform'){
    fmin = pars$fmin
    fmax = pars$fmax
    dtvec_mat          = outer(tvec,tvec,'-')
    tvec_prod_mat      = outer(tvec,tvec,'*')
    ncp_reg_mat       = ((sin(4*pi*fmax*dtvec_mat) - sin(4*pi*fmin*dtvec_mat))*(4*tvec_prod_mat*pi^2 + 1))/(4*dtvec_mat)
    diag(ncp_reg_mat) = pi*(4*pi^2*tvec^2 + 1)*(fmax - fmin) 
    return(sum(sum(Amp^2*ncp_reg_mat)))
  }else{
    stop('unknown choice of prior')
  }
}