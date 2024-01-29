test_that("function conserves inds, moves at least one", {
  Nfine = sample(c(50:100),1)
  Nmeas = sample(c(1,20),1)
  tau      = c(1:Nfine)/Nfine -1/Nfine
  xinds    = sample(c(1:Nfine),Nmeas)
  control  = list(tfun_choice='single-flip')
  expect_equal(length(unique(tfun_svdpower_discrete(xinds,Nfine,control))),Nmeas)
  expect_equal(length(which(abs(xinds-tfun_svdpower_discrete(xinds,Nfine,control))>0)),1)
})
