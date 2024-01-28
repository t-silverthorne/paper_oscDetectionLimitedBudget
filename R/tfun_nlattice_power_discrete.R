#' Requires control$tfun_choice
#' Requires opts$lattice_cstr='sa_lattice' 
tfun_nlattice_power_discrete=function(control,x=NULL,opts){
  if (control$tfun_choice=='default'){
    sa_propfunction(opts,x)
  }
}