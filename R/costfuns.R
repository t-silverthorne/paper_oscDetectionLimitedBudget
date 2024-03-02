#'Cost function for continuous-time power optimisation
#' @param mt vector of measurement times 
#' @param freqs vector of frequencies at which cost function should be evaluated
#' @param Amp amplitude of signal, only relevant if \code{cfuntype=='power'}
#' @param alpha type I error cutoff, only relevant if \code{cfuntype=='power'}
#' @param weight_ncp factor for scaling non-centrality parameter in multi-objective scalarization. 
#' @param regL1 factor for scaling L1-frequency regularisation in multi-objective scalarization. 
#' @param regFder factor for scaling frequency derivative in multi-objective scalarization. 
#' @param gapPenalty factor for scaling adding a penalty for gaps between consecutive measurement times
#' @param leveragePenalty currently depracated, a factor for penalising by statistical leverage
#' @param cfuntype either 'power' or 'ncp' if non-centrality parameter is to be optimised directly
#' 
#' @return minimum of cost function evaluated at \code{freqs}
#' @author Turner Silverthorne
#' @export
costfun_svdpower=function(mt,freqs,
                          Amp=1,alpha=.05,weight_ncp=1,regL1=0,regFder=0,
                          gapPenalty=0,leveragePenalty=0,
                          cfuntype='power'){
  N   = length(mt)
  if (regFder!=0 & cfuntype=='power'){
    stop('Use cfuntype=ncp, frequency regularization for power has not been implemented.')
  }
  if (cfuntype=='ncp'){
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {getMinEig(.,is_symmetric=T)} 
    })
    
  } else if (cfuntype=='power'){
    f0  = qf(p=1-alpha,df1=2,df2=N-3)
    cvals = freqs %>% sapply(function(freq){
    getReducedFIM(mt,list(freq=freq)) %>% {Amp^2*getMinEig(.,is_symmetric=T)} %>% 
                    {1 - pf(q=f0,df1=2,df2=N-3,ncp=.)} 
    })
  }else{
    stop('unknown cfuntype')
  }
  cval = weight_ncp*min(cvals)
  if (regL1>0){
    cval = cval+regL1*mean(cvals)  
  }
  
  if (regFder!=0){
    #dlambda_dfreq = freqs %>% sapply(function(freq){deig_dfreq(mt,freq)}) %>% abs() %>% max() 
    dlambda_dfreq = freqs %>% sapply(function(freq){deig_dfreq(mt,freq)}) %>% {.^2} %>% sum()
    dlambda_dfreq = sqrt(dlambda_dfreq/length(freqs))
    if (is.numeric(regFder)&regFder>0){
      cval = cval-regFder*dlambda_dfreq
    }
    else if(regFder=='2norm'){
      cval = dlambda_dfreq
    }else{
      stop('unrecognised regFder setting')
    }
  }
  
  if (gapPenalty>0){
    cval=cval-gapPenalty*helper_gap_penalty(mt) 
  }
  if (leveragePenalty>0){
    stop('leverage depracated')
  }
  return(cval)
}

#'Cost function for 2-lattice discrete-time power optimisation
#' @description
#' Similar function to [costfun_svdpower] except measurement times can only
#' be collected at union of two uniform grids and these grids are selected from 
#' a discrete underlying uniform grid.
#' 
#' @param N1 number of points in first grid
#' @param dx1 spacing of first grid
#' @param N2 number of points in second grid
#' @param dx2 spacing of second grid
#' @param xshift2 phase shift of second grid
#' @param tau vector of candidate measurement times
#' @param freqs frequencies at which cost function should be evaluated
#'
#' @note This function works in units where \code{dx1} and \code{dx2} are assumed 
#' to refer to indices in the \code{tau} grid.
#' 
#' Remaining arguments are passed to [costfun_svdpower]
#'
#'
#' @return minimum of cost function evaluated at \code{freqs}
#' @author Turner Silverthorne
#' @export
costfun_2lattice_svdpower_discrete = function(N1,dx1,N2,dx2,xshift2,
                                              tau,freqs,...){
xinds=unique(xinds_from_lat1lat2_pars(N1,dx1,N2,dx2,xshift2))
if (any(is.na(xinds))|length(xinds)<N1+N2){
  stop('duplicate or NA meas times present')
}else if(any(xinds>length(tau))){
  stop('edge case')
}else{
  return(costfun_svdpower(tau[xinds],freqs,...))
}
}

#'Cost function for 2-lattice continuous-time power optimisation
#' @description
#' Similar to [costfun_2lattice_svdpower_discrete] except time is no longer restricted
#' to a discrete grid; the lattices can be scaled arbitrarily. 
#' 
#' @param shift1 shift of first grid, must be between 0 and 1
#' @param shift2 shift of second grid, must be between 0 and 1
#' @param scale1 scale factor for first grid
#' @param scale2 scale factor for second grid
#' @param lat1 first uniform grid
#' @param lat2 second uniform grid
#' @param freqs frequencies at which cost function should be evaluated
#' @author Turner Silverthorne
#' @export
costfun_2lattice_svdpower=function(shift1,shift2,scale1,scale2,lat1,lat2,
                                   freqs,...){
  mt=convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2)
  return(costfun_svdpower(mt,freqs,...))
}

#'Cost function for 2-lattice continuous-time power optimisation
#' @description
#' Very similar to [costfun_2lattice_svdpower_discrete] except now, one of the lattices
#' is assumed to be fixed at the roots of unity.
#'
#' This is useful because it removes degeneracy from the problem.
#'
#' @author Turner Silverthorne
#' @export
costfun_auglattice=function(N1,N2,shift2,scale2,freqs,...){
  mt=helper_auglattice_to_state(N1=N1,N2=N2,
                                shift2=shift2,scale2=scale2)
  return(costfun_svdpower(mt,freqs,...))
}

#'Cost function for discrete-time power optimisation
#' @description
#' Similar to [costfun_svdpower] except the points are chosen from a fixed
#' underlying grid, rather than allowed to be arbitrary points in the unit interval.
#'
#' @param xinds the indices of \code{tau} at which measurements should be taken
#' @param tau the underlying fine grid from which measurements should be selected 
#' @param freqs frequencies at which cost function should be evaluated
#' @author Turner Silverthorne
#' @export
costfun_svdpower_discrete = function(xinds,tau,freqs,...){
  return(costfun_svdpower(tau[xinds],freqs,...))
}

