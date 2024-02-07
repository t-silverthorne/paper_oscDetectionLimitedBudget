#' Requires control$tfun_choice
#' Requires opts$lattice_cstr='sa_lattice' 
tfun_nlattice_power_discrete=function(control,x=NULL){
  if (control$tfun_choice=='default'){
    sa_propfunction(control,x)
  }
}