#' Requires control$tfun_choice
#' Requires opts 
tfun_nlattice_svdpower(control,x=NULL,opts){
  if (control$tfun_choice=='default'){
    sa_propfunction(opts,x)
  }
}