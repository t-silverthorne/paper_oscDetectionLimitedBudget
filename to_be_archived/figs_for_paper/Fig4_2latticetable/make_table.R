#WIP  generates data necessary for comparison table
require(dplyr)
require(devtools)
require(ggplot2)
require(data.table)
require(parallel)
require(patchwork)
load_all()

Nmv   = c(16:48)
flist = NULL
for (Nm in Nmv){
  fpin2lat = paste0('figs_for_paper/solutions/solns_pin2lat_',Nm,'.RDS')
  fcvxr    = paste0('figs_for_paper/solutions/solns_cvxr_sparse_',Nm,'.RDS')
  flist=c(fpin2lat,fcvxr,flist) 
}
flist
Nfreq=2^10
freqs = seq(from=1,to=24,length.out=Nfreq)

for (ff in flist){
  dfloc=readRDS(ff)
  if (grepl('cvxr',ff)){
    dfloc$solver='cvxr'
  }
  if (grepl('pin2lat',ff)){
    dfloc$ncp = abs(dfloc$ncp) # consistent sign choice
    maxrep    = dfloc[which.max(dfloc$ncp),]$rep
    dfloc     = dfloc[dfloc$rep==maxrep,] 
    dfloc     = dfloc[, !names(dfloc)%in% c('rep')]
  }
  mtdf=rbind(mtdf,dfloc)
}

res=NULL
for (Nmeas in Nmv){
  mt_unif = c(1:Nmeas)/Nmeas - 1/Nmeas
  mt_cvxr = mtdf[mtdf$solver=='cvxr' & mtdf$Nmeas==Nmeas,]$time
  mt_p2l  = mtdf[mtdf$solver=='pin2lat' & mtdf$Nmeas==Nmeas,]$time
  if (length(mt_cvxr) != Nmeas){
    stop('bad filter')
  }
  if (length(mt_p2l) != Nmeas){
    stop('bad filter')
  }
  
  ncp_unif = costfun_svdpower(mt_unif,freqs,cfuntype = 'ncp')
  ncp_cvxr = costfun_svdpower(mt_cvxr,freqs,cfuntype = 'ncp')
  ncp_p2l  = costfun_svdpower(mt_p2l,freqs,cfuntype = 'ncp')
 
  data.frame(perf_unif  = 100*ncp_unif/ncp_cvxr,
             perf_p2l   = 100*ncp_p2l/ncp_cvxr,
             N1         =
             N2
             ) 
}



