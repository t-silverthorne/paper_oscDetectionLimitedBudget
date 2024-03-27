#'Compute pseudo window of a Lomb Scargle periodogram
#'
#'@description
#' Given a frequency, acrophase, and collection of measurement times, the Lomb-Scargle 
#' periodogram of this signal is referred to as the pseudo-window of the 
#' measurement times
#' 
#' It is useful for estimating how accurately different signals will be detected 
#' by the measurement times.
#' 
#' @param mt vector of measruement times
#' @param freqs frequencies at which to evaluate pseudo window
#' @param acros acrophases at which to evaluate pseudo-window
#' @param type a label for the output, useful for sweeping over multiple collections of measurement times
#'
#' @return a [data.frame] of Lomb-Scargle periodograms for each of the input 
#' frequencies and acrophases
#' 
#' @author Turner Silverthorne
#' @export
lsp_pseudowindow=function(mt,freqs,acros,type='null'){
  Nmeas=length(mt)
  ploc = expand.grid(freq=freqs,acro=acros)
  c(1:dim(ploc)[1]) %>% lapply(function(ind){
    freq=ploc[ind,]$freq
    acro=ploc[ind,]$acro
    pl=lsp(x=cos(2*pi*freq*mt-acro),times=mt,from=fmin,to=fmax,plot=F)
    data.frame(lfreq=pl$scanned,lpower=pl$power,freq=freq,acro=acro,type=type,Nmeas=Nmeas)
  }) %>% rbindlist()
}

