require(devtools)
require(annmatrix)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(annmatrix)
require(lomb)
require(pROC)
require(ggplot2)
sols = readRDS('../hiresSols.RDS')

freqs = c(1)

pars  = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1,2),
                         p_osc = c(0.5),
                         fdr_method=c('none'),
                         type=c('equispaced','threshold','balanced','regu_no_cstr'))