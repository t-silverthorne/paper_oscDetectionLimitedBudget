#' Enumerate partitions of n
#' 
#' @description
#' Use pentagonal number recurrence relation to enumerate partitions of all integers in a range 1 to n
#' 
#' 
#' @author Turner Silverthorne
sa_enumpar=function(N){
  pvec=replicate(N+1,NaN)
  for (n in c(0:N)){
    if (n==0|n==1){
      pvec[n+1]=1
    }else{
      pvec[n+1]=0
      kmin = ceiling(-(sqrt(24*n+1)-1)/6)
      kmax = floor((sqrt(24*n+1)+1)/6)
      S=0
      for (kk in c(kmin:kmax)){
        if (1+n-kk*(3*kk-1)/2 > 0){
          S=S + (-1)^(kk+1)*pvec[1+n-kk*(3*kk-1)/2]
        }
      }
      pvec[n+1]=S
    }
  }
  return(pvec)
}