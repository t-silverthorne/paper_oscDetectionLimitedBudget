#' requires control$tfun_choice
tfun_svdpower_discrete=function(xinds,Nfine,control){
  if (control$tfun_choice=='single-flip'){
    flip_index        = sample(length(xinds),1)
    new_val           = sample(setdiff(c(1:Nfine),xinds),1) 
    xinds[flip_index] = new_val
    #TODO: add discretized Brownian motion
  }else{
    stop('unknown control$tfun_choice')
  }
  return(xinds)
}