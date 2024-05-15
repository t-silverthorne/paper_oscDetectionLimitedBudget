#'Cosinor data generation
#' @description
#' Generates cosinor sample data subject to white noise such as Negative Binomial or Gaussian
#' @param Nmeas Number of measurements 
#' @param Outloop The size of outer loop
#' @param noise_type Either the Negative Binomial or Gaussian (marked as "gaussian" or "negative binomial") which adjusts distribution for white noise
#' @param phase Horizontal translation of the MESOR wave 
#' @param Amp amplitude of the wave 
#' @param period frequency of the peak wave
#' @param parameter vector consisting of amp, phase, and period
#' @return a vector of size Outloop x Nmeas where each row is a schedule of cosinor data and the columns represent each measurement 
#' @author John Anthony Limanto
#' @export

generate_noise = function(noise_type, Outloop, Nmeas){
  if (noise_type == "gaussian") {
    return(matrix(rnorm(Outloop * Nmeas), ncol = Nmeas))
  } else if (noise_type == "negative binomial") {
    return(matrix(rnbinom(outloop * Nmeas, size = 1, prob = 1/2), ncol = Nmeas))
  } else {
    stop("Unsupported noise type")
  }
}

generate_data = function(noise_type, Outloop, Nmeas, param){
  Amp <- param[1]
  phase <- param[2]
  period <- param[3]
  
  mt      = c(1:Nmeas)/Nmeas -1/Nmeas
  cosine_wave <- Amp * cos(2 * pi * mt / period - phase)
  
  Xdat = matrix(rep(cosine_wave, Outloop), nrow = Outloop, byrow = TRUE) +
    generate_noise(noiseType, Outloop, Nmeas)
  
  return(Xdat)
}
