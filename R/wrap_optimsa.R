#' Wrapper for optimising power directly using simulated annealing
#' 
#' @description
#' Wrapper for maximizing worst-case power using [optimization::optim_sa()]. In
#' this function, measurement times are not restricted to a lattice. Worst-case
#' power is interpreted as the minimum power for any frequency between 
#' \code{fmin} and \code{fmax}, and any acrophase in \code{[0,2*pi]}.
#'
#' @param Nmeas measurement budget of experiment
#' @param fmin lowest frequency of interest
#' @param fmax highest frequency of interest
#' @param Nf   number of points to use in freq discretisation
#' @param Nacro number of points to use in acrophase discretisation
#' @param ctrl corresponds to \code{control} parameter of [optimization::optim_sa()]
#' @param Amp amplitude of signal
#' 
#' @return output of simulated annealing
#' 
#' @author Turner Silverthorne
#' @export
wrap_optimsa = function(Nmeas,fmin,fmax,Nf,Nacro,ctrl,Amp=1){
  # unpack parameters
  freqs = seq(from=fmin,to=fmax,length.out=Nf)
  acros = seq(from=0,to=2*pi,length.out=Nphi+1)
  acros = acros[1:Nphi]
  pars  = expand.grid(Amp=Amp,freq=freqs,acro=acros)
 
  # define cost function 
  mt=runif(Nmeas) 
  worst_power = function(mt){
    min(apply(pars,1,function(x){eval_exact_power(t=mt,param=x)}))
  }
  
  # run simulated annealing 
  return(optim_sa(fun=worst_power,
           start=runif(Nmeas),
           lower=replicate(Nmeas,{0}),
           upper=replicate(Nmeas,{1}),
           maximization=TRUE,
           trace = TRUE,
           control=ctrl))
  
}