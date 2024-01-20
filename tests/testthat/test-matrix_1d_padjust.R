library(annmatrix)
library(dplyr)
test_that("dim acts correctly", {
  pvals = runif(10) %>% matrix(nrow=1) %>% as.annmatrix()
  pvals_full = pvals[rep(1,3),]
  
  res1=matrix_1d_padjust(pvals,1,'fdr')
  res2=matrix_1d_padjust(pvals_full,1,'fdr')
  
  # check acting on single row is same as acting on each row
  for (ii in c(1:3)){
    expect_equal(res1[1,],res2[ii,])
  }
  expect_gte(max(res1-pvals),0)
  
  # check dim=2 acts on transpose in expected way
  expect_equal(res2,t(matrix_1d_padjust(t(pvals_full),2,'fdr')))
  
  # check shape is perserved when you use dim=2
  expect_equal(dim(matrix_1d_padjust(pvals_full,1,'fdr')),
               dim(matrix_1d_padjust(pvals_full,2,'fdr')))
})

test_that('row test non-identical data',{
  pvals = runif(30) %>% matrix(nrow=3) %>% as.annmatrix()
  
  res=matrix_1d_padjust(pvals,1,'fdr')
  
  res1=p.adjust(pvals[1,],'fdr')
  res2=p.adjust(pvals[2,],'fdr')
  res3=p.adjust(pvals[3,],'fdr')
  
  expect_equal(res[1,],res1) 
  expect_equal(res[2,],res2) 
  expect_equal(res[3,],res3) 
  
})
