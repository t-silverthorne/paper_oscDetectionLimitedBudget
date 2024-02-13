test_that("returnFder affects output", {
  Nfreq=2^6
  opts=list(Nfine=288,Nfreq=Nfreq,Nmeas=10,fmin=1,fmax=24)
  Amlist=make_quadmats(opts)
  Aderlist=make_quadmats(opts,returnFder=T)

  expect_false(c(1:Nfreq) %>% lapply(function(ii){
    (Amlist[[ii]]==Aderlist[[ii]] )%>% as.numeric() %>% as.logical() %>% any()
    }) %>% unlist() %>% any())
})

test_that("returnFder affects output", {
  Nfreq=1
  df=1e-6
  f0=runif(1)*10
  opts=list(Nfine=100,Nfreq=Nfreq,Nmeas=10,fmin=f0,fmax=f0)
  A  = make_quadmats(opts) %>% {.[[1]]}
  dAdf=make_quadmats(opts,returnFder=T) %>% {.[[1]]}
  
  opts$fmin=f0+df
  opts$fmax=f0+df
  Ah = make_quadmats(opts) %>% {.[[1]]}

  dAdf_num = (Ah-A)/df
  
  expect_lt(max(abs(dAdf-dAdf_num))/max(abs(dAdf_num)),10*df)
})