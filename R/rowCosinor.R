#' Cosinor model
#'
#' Function to fit a cosinor model for each row of the matrix.
#'
#' Fits a cosinor model without any additional vocatiates.
#'
#' @param theData input matrix
#' @param zts time values for each column of theData
#' @param per period of oscillation to be fitted (default is 24)
#'
#' @return a data.frame containing various results of the fit including
#' p-values, R-squared values, acrophases and amplitudes for each row.
#'
#' @author Karolis Koncevicius
#' @export
rowCosinor <- function(theData, zts, per=24) {
  Y <- as.matrix(theData)

  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)

  betas <- solve(t(X) %*% X) %*% t(X) %*% t(Y)

  phases     <- atan2(betas[2,], betas[3,])
  acrophases <- (((per/2) / pi) * (phases)) %% per
  amplitudes <- sqrt(betas[2,]*betas[2,] + betas[3,]*betas[3,])

  fits <- t(X %*% betas)

  SStot <- rowSums((Y - rowMeans(Y))^2)
  SSres <- rowSums((fits-Y)^2)
  Rsqs  <- 1 - (SSres/SStot)

  SSmod <- SStot - SSres
  DFres <- ncol(theData) - 3
  DFmod <- 2
  MSres <- SSres / DFres
  MSmod <- SSmod / DFmod
  Fstatistic <- MSmod / MSres

  pval <- pf(Fstatistic, DFmod, DFres, lower.tail=FALSE)

  data.frame(acrophase=acrophases, amplitude=amplitudes, mesor=betas[1,],
             rsq=Rsqs, statistic=Fstatistic, pvalue=pval
             )
}

