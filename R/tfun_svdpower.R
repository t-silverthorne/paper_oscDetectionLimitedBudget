#' requires control$tfun_choice
#'          control$tfun_mean
#'          control$tfun_sd 
tfun_svdpower=function(mt,control){
  if(control$tfun_choice=='brownian'){
    mt=mt+rnorm(length(mt),mean=control$tfun_mean,sd = control$tfun_sd) 
  }else{
    stop('unknown control$tfun_choice')
  }
  return(mt)
}