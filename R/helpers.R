#'Convert augmented lattice to a vector of measurement times
#' 
#' @param N1 number of points in first lattice
#' @param N2 number of points in second lattice
#' @param shift2 offset of second lattice
#' @param scale2 scale factor for second lattice
#' @return vector of measurement times corresponding to the input 2-lattice parameterisation 
#' 
#' @author Turner Silverthorne
#' @export
helper_auglattice_to_state=function(N1,N2,shift2,scale2){
  lat1=c(1:N1)/N1-1/N1
  lat2=c(1:N2)/N2-1/N1
  lat2=shift2+scale2*lat2
  mt=c(lat1,lat2)
  return(mt)
}

#' Helper function converts parameterization of two lattices into single lattice
#' 
#' @param shift1 shift of first lattice 
#' @param shift2 shift of second lattice 
#' @param scale1 scale of first lattice
#' @param scale2 scale of second lattice
#' @param lat1 vector of measurement times in first lattice
#' @param lat2 vector of measurement times in second lattice
#' 
#' @return vector of measurement times corresponding to the input 2-lattice parameterisation 
#' 
#' @author Turner Silverthorne
#' @export
convert_2lattice_to_state=function(shift1,shift2,scale1,scale2,lat1,lat2){
  x1 = shift1+scale1*lat1
  x2 = shift2+scale2*lat2
  # keep measurement times inside study, necessary for L-BFGS-B solver but not 
  # simulated annealing (in simulated annealing, times stay inside study because
  # of design of transition function) 
  x  = c(x1,x2) 
  x  = x %% 1
  return(x)
}

#' Helper function for lattice-based simulated annealing
#' 
#' @param scale current scale of lattice 
#' @param shift current state of lattice 
#' @param lat vector of measurement times in current lattice 
#' @param tscale width parameter for uniform random points 
#' 
#' @return updates the shift of the lattice by sampling a uniform distribution 
#' of width \code{tscale}, rejection sampling is then used to make sure that the
#' lattice with the updated shift parameter is still contained in the unit interval.
#' 
#' @author Turner Silverthorne
#' @export
helper_unif_update_shift=function(scale,shift,lat,tscale){
  if (any(scale*lat + shift <0)  | any(scale*lat + shift >1)){
    stop('invalid initial state')
  }else{
    # TODO: replace with simpler parameterization 
    lmin=min(lat)
    lmax=max(lat) 
   
    sh_min = -scale*lmin
    sh_max = 1-scale*lmax
     
    eps_n = min(tscale/2,shift-sh_min)
    eps_p = min(tscale/2,sh_max-shift)
   
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    }
    shift_new = runif(1,shift-eps_n,shift+eps_p)
    return(shift_new)
  }
}

#' Discretised version of shift update
#' 
#' @seealso [helper_unif_update_shift]
#' @author Turner Silverthorne
#' @export
helper_unif_update_shift_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (scale*(Npts-1)+shift > Nfine | shift<0 | scale < 0){
    stop('invalid initial state')
  }else{
    sh_min = 0
    sh_max = Nfine-1 -(Npts-1)*scale
    
    xL = max(sh_min,shift-floor(tscale/2))
    xR = min(sh_max,shift+floor(tscale/2))
    
    shift_new = sample(c(xL:xR),1)
    return(shift_new)
  }
}


#' Helper function for lattice-based simulated annealing
#' 
#' @param scale current scale of lattice 
#' @param shift current state of lattice 
#' @param lat vector of measurement times in current lattice 
#' @param tscale width parameter for uniform random points 
#' 
#' @return updates the scale of the lattice by sampling a uniform distribution 
#' of width \code{tscale}, rejection sampling is then used to make sure that the
#' lattice with the updated scale parameter is still contained in the unit interval.
#' 
#' @author Turner Silverthorne
#' @export
helper_unif_update_scale=function(scale,shift,lat,tscale){
  if (any(scale*lat + shift <0)  | any(scale*lat + shift >1)){
    stop('invalid initial state')
  }else{
    # TODO: replace with simpler parameterization 
    lmin=min(lat)
    lmax=max(lat) 
   
    sc_min = 0
    sc_max = (1-shift)/lmax
     
    eps_n = min(tscale/2,scale-sc_min)
    eps_p = min(tscale/2,sc_max-scale)
    if (eps_n <0 | eps_p < 0){
      stop('negative window')
    } 
    scale_new = runif(1,scale-eps_n,scale+eps_p)
    return(scale_new)
  }
}


#' Discretised version of scale update
#' 
#' @seealso [helper_unif_update_scale]
#' @author Turner Silverthorne
#' @export
helper_unif_update_scale_discrete=function(scale,shift,Npts,tscale,Nfine){
  if (tscale<2){
    stop('for discrete problem, tscale must be >= 2')
  }
  if (scale*(Npts-1)+shift > Nfine | shift<0 ){
    stop('invalid initial state')
  }else{
    sc_min = 1 
    sc_max = floor((Nfine-shift-1)/(Npts-1)) #TODO: verify you don't need floor here
    
    xL = max(sc_min,scale-floor(tscale/2))
    xR = min(sc_max,scale+floor(tscale/2))
    
    scale_new = sample(c(xL:xR),1)
  }
  return(scale_new)
}

#'For penalizing cost function based on consecutive measurement gap 
#' 
#' @param mt vector of measurement times in unit interval \code{[0,1]}
#' 
#' @return the largest gap between consecutive measurements.
#' 
#' @note This includes the gap between the final measurement and first measurement
#' \code{1-(mt[length(mt)]-mt[1])}
#' 
#' @author Turner Silverthorne
#' @export
helper_gap_penalty=function(mt){
  mt=sort(unique(mt))
  if (any(mt>1) | any(mt<0)){
    stop('measurement times must be in [0,1]')
  }
  d1 = 1-(mt[length(mt)]-mt[1])
  return(max(c(diff(mt),d1)))
}

#' Helper function for discrete two-lattice optimization
#' @description
#' Given a fine-uniform grid and a parameterization of a 2-lattice contained
#' in this grid, this function returns the grid indicest that correspond to this parameterization.
#' 
#' @param N1 number of points in first lattice
#' @param dx1 spacing of first lattice 
#' @param N2 number of points in second lattice
#' @param dx2 spacing of second lattice 
#' @param xshift2 offset of second lattice
#' 
#' @return a vector of indices corresponding to measurement times
#' 
#' @note all parameters should be integers, as they are parameterising indicies
#' 
#' @author Turner Silverthorne
#' @export
xinds_from_lat1lat2_pars=function(N1,dx1,N2,dx2,xshift2){
  lat1  = 1+dx1*c(0:(N1-1)) 
  lat2  = 1+xshift2 + dx2*c(0:(N2-1))
  xinds = c(lat1,lat2)
  return(xinds)
}

#'Sweep over different partiions of a measurement budget and return optimal
#'
#' @description
#' Given an \code{Nmeas} budget, split \code{Nmeas=N1+N2} for various values of 
#' \code{N1} and \code{N2} and perform optimization in this subspaces
#' 
#' @param optim_routine the optimization routine that should be run
#' @param Nvals which values of N1 should be considered
#' @param Nmeas total measurement budget
#' @param cts_flag is your optimization routine using continuous time?
#' @param dvar0 initial state of optimizer
#' 
#' @return a list containing the optimal solution, total computation time, and the
#' values of the objective function at all other Nvals for comparison.
#' @author Turner Silverthorne
#' @export
wrapper_sweep_lattice=function(optim_routine,Nvals,Nmeas,cts_flag,dvar0,...){
  fval_best = Inf
  res_best  = NULL
  N1_best   = NULL
  time_tot  = 0
  all_fvals = NULL 
  for (n1 in Nvals){
    n2=Nmeas-n1
    if (cts_flag){
      dvar0[['lat1']]= c(1:n1)/n1 - 1/n1
      dvar0[['lat2']]= c(1:n2)/n2 - 1/n2
      res_loc = optim_routine(dvar0,...)
    }else{
      dvar0$N1=n1
      dvar0$N2=n2
      res_loc = optim_routine(dvar0,...)
    }
    time_tot=time_tot+as.numeric(res_loc$runtime)
    all_fvals=rbind(all_fvals,data.frame(fval=res_loc$fvalue,N1=n1))
    if (res_loc$fvalue<fval_best){
      res_best  = res_loc
      N1_best   = n1
      fval_best = res_loc$fvalue
    }
  }
  return(list(res_best=res_best,N1_best=N1_best,
              time_tot=time_tot,
              all_fvals=all_fvals))
}

