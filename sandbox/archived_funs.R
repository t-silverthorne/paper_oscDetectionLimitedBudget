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


#'Generate decision variable for convex power optimization 
#' @description
#' generates a binary vector which represents the weight function used by CVXR
#' for power optimization.
#' 
#' @param opts$Nfine the number of points in the underlying fine grid from which 
#' measurement times are to be chosen.
#'   
#' 
#'   
#' @return \code{x} binary [CVXR::Variable] to be used as the decision variable in CVXR programming.
#' @author Turner Silverthorne
#' @export
make_variable=function(opts){
  x=Variable(opts$Nfine,boolean=T)
  return(x)
}

#' Generate matrices for quadratic forms used in convex power optimization
#' 
#' @param opts$Nfine number of points in the fine time grid
#' @param opts$Nfreq number of frequencies in frequency grid
#' @param opts$Nmeas measurement budget
#' @param opts$fmin minimum frequency to be included
#' @param opts$fmax maximum frequency to be included
#' 
#' @return a list of matrices of length \code{opts$Nfreq}. Each entry in the list 
#' corresponds to a quadratic form at a frequency between \code{opts$fmin} and 
#' \code{opts$fmax}.
#' 
#' @author Turner Silverthorne
#' @export
make_quadmats = function(opts,...){
  Nfine   = opts$Nfine
  Nfreq   = opts$Nfreq
  Nmeas   = opts$Nmeas
  fmin    = opts$fmin
  fmax    = opts$fmax
  
  tau   = (0:Nfine)/Nfine       
  tau   = tau[1:(length(tau)-1)] 
  fvec  = seq(from=fmin,to=fmax,length.out=Nfreq)
  
  Amlist = list() 
  for (ii in c(1:Nfreq)){
    freq  = fvec[ii]
    cvec  = matrix(cos(2*pi*freq*tau),nrow=Nfine)
    svec  = matrix(sin(2*pi*freq*tau),nrow=Nfine)
    
    a11   = cvec*cvec
    a22   = svec*svec
    a12   = cvec*svec
    
    Amat  = a11%*%t(a11) + a22%*%t(a22) +4*a12%*%t(a12)-a11%*%t(a22) - a22%*%t(a11)
    
    Amlist[[ii]]=Amat 
    
  }

  return(Amlist)  
}

#' Generate a convex programming problem equivalent to power optimization
#' 
#' @description
#' generates a disciplined convex programming problem equivalent to maximizing 
#' the statistical power of a cosinor-based hypothesis test. The problem is created 
#' as a [CVXR::Problem()] which can then be solved using [CVXR] together with the 
#' [GUROBI] backend. The [GUROBI] backend is necessary because there are integer
#' constraints on the problem.
#' 
#' @param x a [CVXR::Variable()] representing the weight function for measurement times, see [make_variable()]
#' @param Aquad a list of symmetric matrices generated by [make_quadmats]
#' @param opts$lattice_cstr must be equal to \code{'none'} for function to evaluate 
#' @param opts$optim_method must be equal to \code{'cvxr'} for function to evaluate 
#' 
#' @return \code{prob} a disciplined convex programming problem [CVXR::Problem] 
#' which can be solved using the [CVXR] library. 
#' 
#' @author Turner Silverthorne
#' @export
make_problem=function(x,Aquad,opts,regL1=0,regFder=0,...){
  # convert to CVXR datatype
  for (ii in c(1:length(Aquad))){
    Aquad[[ii]]=Constant(Aquad[[ii]])
  }
  
  Nm   = Constant(opts$Nmeas)
  csts = list( sum(x) == Nm)
  big_str = paste(lapply(1:length(Aquad) %>% as.list(),
         function(ind){
           paste0('quad_form(x,Aquad[[',ind,']])')
         }),collapse=',')
 
  str_prefix='prob=Problem(Minimize('
  str_suffix=paste0('max_elemwise(',big_str,')),csts)')
  
  if(regL1>0){
    regL1 = Constant(regL1)  
    str_prefix=paste0(str_prefix,'regL1*quad_form(x,Amean)+')
  }
  if(regFder>0){
    stop('freq regularization not yet implemented for CVXR')
  }
  
  strp=paste0(str_prefix,str_suffix)
  eval(parse(text=strp))
  return(prob)
}

make_problem_entries=function(x,Aquad,opts){
  Nm   = Constant(opts$Nmeas)
  csts = list( sum(x) == Nm)
  for (ii in c(1:length(Aquad))){
    if (ii == 1){
      obj = quad_form(x,Constant(Aquad[[ii]]))
    }
    obj = vstack(obj,quad_form(x,Constant(Aquad[[ii]])))
  }
  
  prob = Problem(Minimize(max_entries(obj)),csts)
}


make_problem_hvolume=function(x,Aquad,opts,weight){
  Nm   = Constant(opts$Nmeas)
  csts = list( sum(x) == Nm)
 
  Amean_w= 0
  for (ii in c(1:length(Aquad))){
    Amean_w= Amean_w + Aquad[[ii]]*weight[ii]
  }
  Amean_w=Amean_w/length(Aquad)
  Amean_w=Constant(Amean_w)
 
  
  prob=Problem(Minimize(quad_form(x,Amean_w)),csts)
  return(prob)
}

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
solve_cvxr_spt= function(control,Threads,use_spt=F,drts=NULL,pads=NULL){
  Aquad = make_quadmats(control)
  x     = Variable(control$Nfine,boolean=T)
  
  #  ensure support distributed throughout domain
  csts  = list(sum(x)==control$Nmeas)
  
  if (use_spt){
    rt_inds=seq(1,control$Nfine,drts)
    for (ii in c(1:length(rt_inds))){
      cst_spt = sum(x[c(rt_inds[ii]:(rt_inds[ii]+pads))])>=1
      csts = append(csts,cst_spt)
    }
  }
  
  # convert to CVXR datatype
  for (ii in c(1:length(Aquad))){
    Aquad[[ii]]=Constant(Aquad[[ii]])
  }
  
  str_prefix='prob=Problem(Minimize('
  big_str = paste(lapply(1:length(Aquad) %>% as.list(),
                         function(ind){
                           paste0('quad_form(x,Aquad[[',ind,']])')
                         }),collapse=',')
  
  str_suffix=paste0('max_elemwise(',big_str,')),csts)')
  
  strp=paste0(str_prefix,str_suffix)
  eval(parse(text=strp))
  
  xout=CVXR::solve(prob,verbose=control$cvxr_verbose,num_iter=control$maxit,solver="GUROBI",
                   TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,MIPGap=control$MIPGap,
                   Presolve=control$PreSolve,MIPFocus=control$MIPFocus,
                   Threads=Threads,NodefileStart=control$NodefileStart)
  #tau = c(1:control$Nfine)/control$Nfine - 1/control$Nfine
  #mtvalue    = tau[as.logical(xout[[1]]>1-1e-6)] # TODO: better way of catching this 
  
  return(xout)
}
#'Solve power optimization using CVXR
#' 
#' @description
#' Performs DCP if \code{control$prob_formulation = 'max_elemwise'} or \code{'max_entries'}.
#' Performs hypervolume scalarization if \code{control$prob_formulation == 'hvolume'}.
#' 
#' @param control$Nmeas number of measurements to be collected
#' @param control$Nfine number of points in fine grid from which measurement times are to be chosen
#' @param control$Nfreq number of points to use in frequency discretisation of objective
#' @param control$fmin minimum candidate frequency
#' @param control$fmax maximum candidate frequency
#' @param control$PreSolve recommended \code{PreSolve=2}, see [CVXR::solve()] for details
#' @param control$MIPFocus recommended \code{PreSolve=3}, see [CVXR::solve()] for details
#' @param control$cvxr_verbose boolean for how much output is desired from CVXR
#' @param control$num_iter max number of CVXR iterations
#' @param control$MIPGapAbs stopping criterion for CVXR solver
#' @param control$time_limit max time allowed for computation
#' @param control$prob_formualtion either max_elemwise, max_entries, or hvolume, see description
#' @param Threads how many threads should CVXR use 
#'  
#' @return the result of calling [CVXR::solve] with specified user controls
#' 
#' @note
#' the \code{...} argument is only relevant for computing the cost function value
#' at the end of this method, it is not relevant to the optimization stage.
#' 
#' @note 
#' \code{time_limit} does not count the time taken to complete the pre-solve.
#' 
#' @note 
#' if iterations exceed \code{num_iter}, it is difficultto access output of the computation.
#' For this reason, we suggest setting \code{num_iter} to an unreasonably high 
#' value and using \code{time_limit} instead to determine the cut-off time for computation.
#' 
#' @seealso [opt_osc_power()] which contains a similar wrapper for CVXR
#' @author Turner Silverthorne
#' @export
solve_cvxr=function(control,Threads,...){
  start_time=Sys.time()
  weight=NULL
  Aquad             = make_quadmats(control)
  x                 = make_variable(control)
  if (control$prob_formulation=='max_elemwise'){
    prob = make_problem(x,Aquad,control)
  }else if (control$prob_formulation=='max_entries'){
    prob = make_problem_entries(x,Aquad,control) 
  }else if (control$prob_formulation=='hvolume'){
    weight=rnorm(length(Aquad)) # sample positive orthant of unit sphere
    weight=abs(weight/sqrt(sum(weight^2))) #TODO: verify
    prob = make_problem_hvolume(x,Aquad,control,weight)
  }else if(control$prob_formulation=='L1'){
    weight=rep(1,length(Aquad))/length(Aquad)
    prob = make_problem_hvolume(x,Aquad,control,weight)    
  }else{
    stop('unrecognised problem formulation')
  }
  xout = CVXR::solve(prob,verbose=control$cvxr_verbose,num_iter=control$maxit,solver="GUROBI",
                       TimeLimit=control$time_limit,MIPGapAbs=control$MIPGapAbs,
                       Presolve=control$PreSolve,MIPFocus=control$MIPFocus,
                     Threads=Threads,NodefileStart=control$NodefileStart)
  
  tau = c(1:control$Nfine)/control$Nfine - 1/control$Nfine
  mtvalue    = tau[as.logical(xout[[1]]>1-1e-6)] # TODO: better way of catching this 
  
  freqs      = seq(from=control$fmin,to=control$fmax,length.out=control$Nfreq)
  fvalue     = costfun_svdpower(mt=mtvalue,freqs=freqs,...) 
  xindsvalue = xout[[1]]
  
  tstamp   = lubridate::now() %>% toString() %>% str_replace(' ','___')
  
  end_time =Sys.time()
  
  res_full = list(fvalue     = fvalue,
       mtvalue    = mtvalue,
       xindsvalue = xindsvalue,
       timestamp  = tstamp,
       runtime    = end_time-start_time,
       optim_raw  = xout,
       weight     = weight)
}