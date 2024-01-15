test_that("function evaluates", {
  fmin  = 1
  fmax  = 20
  Nf    = 2^2
  Nacro = 2^5
  Nmeas = 12
  nlim=2e2
  Amp=1.6
  ctrl = list(nlimit=nlim,
              r=.95,
              dyn_rf=T,
              ac_acc=.1,
              stopac=nlim
              )
  uu=wrap_optimsa(Nmeas,fmin,fmax,Nf,Nacro,ctrl,Amp)
  plot(uu)
})

test_that("function evaluates", {
  fmin  = 1
  fmax  = 20
  Nf    = 2^2
  Nacro = 2^5
  Nmeas = 12
  nlim=2e2
  Amp=1.6
  vf_cust =function(para_0){
    ret_var_func=para_i+.01*rnorm(length(para_i))
    return(ret_var_func)}
  var_func <- function(para_0, fun_length, rf) {
    ret_var_func <- para_0 + runif(fun_length, 1e-06, 
                                   rf) * ((rbinom(fun_length, 1, 0.5) * -2) + 1)
    return(ret_var_func)}
  ctrl = list(nlimit=nlim,
              r=.1,
              dyn_rf=T,
              ac_acc=.1,
              stopac=nlim,
              vf=var_func
              )
  uu=wrap_optimsa(Nmeas,fmin,fmax,Nf,Nacro,ctrl,Amp)
  plot(uu)
})
