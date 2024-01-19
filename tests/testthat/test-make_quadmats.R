library(CVXR)

test_that('Aquad non-negative def L1',{
  opts = make_default_opts(prob_size='medium',
                           solver_type='simulanneal')
  opts$costfun_type = 'L1'
  Amat = make_quadmats(opts)
 
  evals=eigen(Amat,symmetric=T) %>% {.$values}
  min(evals)
  expect_gte(min(evals),-1e-13)
})


test_that('Aquad non-negative def Linfty',{
  opts = make_default_opts(prob_size='medium',
                           solver_type='simulanneal')
  opts$costfun_type = 'Linfty'
  Amat = make_quadmats(opts)

  length(Amat) 
  for (ii in c(1:length(Amat))){
    Aloc = Amat[[ii]]
    evals=eigen(Aloc,symmetric=T) %>% {.$values}
    min(evals)
    expect_gte(min(evals),-1e-13)
  }
})

test_that('Aquad symmetric',{
  
  opts = make_default_opts(prob_size='medium',
                           solver_type='simulanneal')
  opts$costfun_type = 'Linfty'
  Amat = make_quadmats(opts)

  length(Amat) 
  for (ii in c(1:length(Amat))){
    Aloc = Amat[[ii]]
    expect_lte(max(abs(Aloc-t(Aloc))),1e-15)
  }
})