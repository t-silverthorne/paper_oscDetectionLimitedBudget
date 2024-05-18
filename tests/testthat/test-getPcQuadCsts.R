test_that("regularization only matches", {
  fmin=1.15
  fmax=24.1
  
  QC=getPcQuadCsts(fmin,fmax,100,w_reg=1,w_wc=0,3) 
  A=matrix(c(
    0.072099551399886,   0.001050943903833,   0.000356646193611,
    0.001050943903833,   0.388363573538944,   0.010270855529227,
    0.000356646193611,   0.010270855529227,   1.337155639956118
  ),nrow=3)
  expect_equal(QC$Qc[c(1:3),c(1:3)],1000*A,tolerance =1e-6)
})