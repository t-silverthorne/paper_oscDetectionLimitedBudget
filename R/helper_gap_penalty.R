helper_gap_penalty=function(mt){
  mt=sort(unique(mt))
  if (any(mt>1) | any(mt<0)){
    stop('measurement times must be in [0,1]')
  }
  d1 = 1-(mt[length(mt)]-mt[1])
  return(max(c(diff(mt),d1)))
}