#' Helper function for sa_propfunction
#'
#' @author Turner Silverthorne

sa_test_for_overlap=function(xnew,x,inds,opts){
  if(opts$enforce_overlap =='ignore'){
    overlap_cond_met=T
  }else if(opts$enforce_overlap =='use-reject'){
    if (length(x)==0){
        ytmp = xnew
    }else{
      ytmp = append(x[!(c(1:length(x)) %in% inds)],xnew) 
    }
    meas_sched = Reduce('+',ytmp)
    overlap_cond_met = as.logical(prod(meas_sched%in%c(0,1)))
  }else{
    stop('opts$enforce_overlap not recognized, please specify valid option')
  }
  return(overlap_cond_met)
}