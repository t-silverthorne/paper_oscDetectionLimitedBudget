#' requires control$tfun_choice
#'          control$mean
#'          control$sd 
tfun_svdpower=function(mt,control){
  if(control$tfun_choice=='brownian'){
    mt=mt+rnorm(length(mt),mean=control$mean,sd = control$sd) 
  }else{
    stop('unknown control$tfun_choice')
  }
  return(mt)
}