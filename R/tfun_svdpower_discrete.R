#' requires control$tfun_choice
tfun_svdpower_discrete=function(xinds,control){
  if (control$tfun_choice=='single-flip'){
    turn_on         = sample(which(xinds==0),1)
    turn_off        = sample(which(xinds==1),1)
    xinds[turn_on]  = 1
    xinds[turn_off] = 0
  }else{
    stop('unknown control$tfun_choice')
  }
  return(xinds)
}