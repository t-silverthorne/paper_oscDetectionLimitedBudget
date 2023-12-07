#' Cost function for simulated annealing power
#' 
#' 
#' @author Turner Silverthorne

cfun_sim_anneal_pwr=function(x,Amat,opts){
  if (opts$costfun_type == 'Linfty'){
    if ('list' %in% class(Aquad)){
      return(max(unlist(lapply(Amat,function(A){t(x)%*%A%*%x}))))
    }else{
      stop("For Linfty cost fun, Aquad should be a list of matrices not a single matrix")
    }
  }else if(opts$costfun_type=='L1'){
    if ('list' %in% class(Aquad)){
      stop("For L1 cost fun, Aquad should be a single matrix not a list of matrices")
    }else{
      return(t(x)%*%Amat%*%x)
    } 
  }else{
    stop('Cost function name not recognized, check opts$costfun_type')
  }
}