test_that("error unrec lprop", {
  opts=list(tfun_choice ='foo')
  expect_error(sa_sample_inds(10,opts),'Unrec')
})

test_that("extremal test linear", {
  m    = 8
  opts = list(tfun_choice='linear')
  dat  = replicate(1e3,{length(sa_sample_inds(m,opts))})
  expect_gte((dat==1 )%>% sum(),(dat==m )%>% sum())
})

test_that("extremal test exp", {
  m    = 8
  opts = list(tfun_choice='exp')
  dat  = replicate(1e3,{length(sa_sample_inds(m,opts))})
  expect_gte((dat==1 )%>% sum(),(dat==m )%>% sum())
})

test_that("exp linear compare", {
  m    = 8
  opts = list(tfun_choice='linear')
  dat1 = replicate(1e3,{length(sa_sample_inds(m,opts))})
  
  opts = list(tfun_choice='exp')
  dat2 = replicate(1e3,{length(sa_sample_inds(m,opts))})
  expect_lte((dat1==1 )%>% sum(),(dat2==1 )%>% sum())
})

#TODO: more precise check of sampling frequencies




