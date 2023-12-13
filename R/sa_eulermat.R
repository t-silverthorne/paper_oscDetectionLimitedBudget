#' Generate the probability matrices used for generating random partitions
#' 
#' @description
#' Helper function for [sa_randpar()], generates the probability matrices used in efficiently
#' sampling random partitions of a given integer using a combinatorial identity of Euler.
#' 
#' @author Turner Silverthorne
sa_eulermat = function(m){
  Pmat = list()
  pvec = sa_enumpar(m)
  for (jj in c(1:m)){
    for (dd in c(1:m)){
      if (jj*dd <= m){
        Pmat[[length(Pmat)+1]]= c(dd,jj,dd*pvec[1 + m - jj*dd]/m/pvec[1+m])
      }
    }
  }
  return(Pmat)
}
