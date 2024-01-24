tfun_svdpower=function(mt,control){
  if(control$tfun_choice=='brownian'){
    mt=mt+rnorm(length(mt),mean=control$mean,sd = control$sd) 
  }
  return(mt)
}