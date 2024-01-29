test_that("points stay inside lattice", {
  control=list(tfun_choice='default')
  control$mean_scale1=0
  control$mean_scale2=0
  control$mean_shift2=0
  control$sd_scale1=.1
  control$sd_scale2=.1
  control$sd_shift2=.1
  shift1=runif(1)
  shift2=runif(1)
  scale1=runif(1)
  scale2=runif(1)
  N1=sample(c(8:12),1)
  N2=sample(c(8:12),1)
  lat1= c(1:N1)/N1 -1/N1
  lat2= c(1:N2)/N2 -1/N2
  for (ii in c(1:100)){
    L=tfun_2lattice_svdpower(shift1,shift2,scale1,scale2,lat1,lat2,control)
    scale1=L$scale1
    scale2=L$scale2
    shift2=L$shift2
  }
  mt=convert_2lattice_to_state(shift1,shift2,scale1,scale2,lat1,lat2)
  expect_false(any(mt>1))
  expect_false(any(mt<0))
})
