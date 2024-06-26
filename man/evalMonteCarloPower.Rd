% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evalMonteCarloPower.R
\name{evalMonteCarloPower}
\alias{evalMonteCarloPower}
\title{Use Monte Carlo method to estimate power}
\usage{
evalMonteCarloPower(
  tvec,
  param,
  Nmc,
  alpha = 0.05,
  method = "Ftest",
  Nperm = 100
)
}
\arguments{
\item{tvec}{vector of measurement times}

\item{Nmc}{number of Monte Carlo samples}

\item{alpha}{type I error, default value \code{alpha=.05}}

\item{param$Amp}{amplitude of signal}

\item{param$freq}{frequency fo signal}

\item{param$acro}{phase of signal in radians}
}
\value{
Monte Carlo estimate of power
}
\description{
Performs harmonic regression on a simulated dataset of independent samples
and returns the portion of samples that have statistically significant p-values.
Useful for comparison with the exact expression for statistical power,
see \link{eval_exact_power}. Relies on \link[matrixTests:cosinor]{matrixTests::row_cosinor}.
}
\author{
Turner Silverthorne
}
