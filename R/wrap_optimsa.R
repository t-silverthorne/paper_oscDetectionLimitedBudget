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
#' @param Amp amplitude of signal
#' 
#' @return output of simulated annealing
#' 
#' @author Turner Silverthorne
#' @export
wrap_optimsa = function(Nmeas,fmin,fmax,Nf,Nacro,Amp=1){
  # unpack parameters
  freqs = seq(from=fmin,to=fmax,length.out=Nf)
  acros = seq(from=0,to=2*pi,length.out=Nacro+1)
  acros = acros[1:Nacro]
  pars  = expand.grid(Amp=Amp,freq=freqs,acro=acros)
 
  # define cost function 
  mt=runif(Nmeas) 
  worst_power = function(mt){
    -min(apply(pars,1,function(x){eval_exact_power(t=mt,param=x)}))
  }
  
  #TODO: make these inputs
  verbose      = T 
  epochs       = F
  num_iter     = 3e4 
  cooling_rate = .999 
  Tinit        = 1000
  Tnow         = Tinit
  rands        = runif(num_iter)
  Sout         = replicate(num_iter,{NaN}) 
  Tout         = replicate(num_iter,{NaN}) 
  ii=1
  x = runif(Nmeas)
  Sx=worst_power(x)
  while (ii<num_iter+1){
    if (ii %% ceiling(num_iter/20) == 0){
      if (verbose){
        print(paste0('Completed: ',
                     toString(100*ii/num_iter),
                     ' perc of SA run. Temp: ',
                     toString(Tnow),
                     ' fval: ',toString(Sx)))
      }
      if (epochs){
        Tnow  = Tinit*.65
        Tinit = Tnow
      }
    }
    y = (x + .01*rnorm(length(x)))%%1 # brownian motion on circle 
    
    Sy=worst_power(y)      # eval cost fun of candidate state
    alpha = min(exp(-(Sy-Sx)/Tnow),1)         # acceptance prob
    
    Tnow = Tnow*cooling_rate
    if (rands[ii]<alpha){
      x  = y
      Sx = Sy
    }
    Sout[ii]=Sx
    Tout[ii]=Tnow
    ii=ii+1
  }
 
  return(list(x=x,Sout=Sout,Tout=Tout))
}