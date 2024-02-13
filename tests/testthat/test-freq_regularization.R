test_that("multiplication works", {
  mt=runif(16)
  
  # cost function value without any regularization
  c1=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),cfuntype = 'ncp')
  c1b=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=0,cfuntype = 'ncp')
  c1c=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regFder=0,cfuntype = 'ncp')
  c1d=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=0,regFder=0,cfuntype = 'ncp')
  expect_equal(c1,c1b)
  expect_equal(c1b,c1c)
  expect_equal(c1c,c1d)
 
  # add mean freq regularization 
  c2=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=1,cfuntype = 'ncp')
  expect_gt(c2,c1)
 
  # add frequency derivative regularization 
  c3=costfun_svdpower(mt,freqs=seq(from=1,to=24,length.out=2^6),regL1=1,regFder=1,cfuntype = 'ncp')
  expect_lt(c3,c2)
  
  expect_error(costfun_svdpower(mt,freqs=c(1,2,3),regFder=10),'Use cfuntype')
})
