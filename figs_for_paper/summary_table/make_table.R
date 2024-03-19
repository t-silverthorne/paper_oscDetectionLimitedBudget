require(kableExtra)
require(dplyr)
require(ggplot2)
require(annmatrix)
require(devtools)
require(data.table)
require(parallel)
load_all()
res=NULL
Nfreq = 2^10
freqs = seq(1,24,length.out=Nfreq)
sol_dir = 'transfer_fold/'
extract_cvxr = function(fname,sol_name){
  dat       = readRDS(paste0(sol_dir,fname))
  xinds     = dat[[1]] 
  Nfine     = length(xinds)
  Nmeas     = sum(xinds)
  tau       = c(1:Nfine)/Nfine -1/Nfine
  mt        = tau[as.logical(xinds)]
  ncp       = costfun_svdpower(mt,freqs,Amp=NaN,cfuntype = 'ncp')
  stat_code = attr(dat$status,'gurobi_status_code')
  stat_desc = attr(dat$status,'gurobi_status_desc')
  data.frame(Nmeas=Nmeas,ncp=ncp,stat_code=stat_code,stat_desc=stat_desc,
             solver=sol_name)
}
extract_cvxr_files = function(fname_list,sol_name){
 return( c(1:length(fname_list)) %>% mclapply(function(ii){
  fname = fname_list[ii]
  extract_cvxr(fname,sol_name)
}) %>% rbindlist() %>% data.frame())
}


# get CVXR solutions
cvxr_files = list.files(sol_dir,pattern = 'solns_cvxr*')
res_cvxr   = extract_cvxr_files(cvxr_files,'cvxr') 


# get CVXR supp solutions
cvxr_supp_files_6 = list.files(sol_dir,pattern = 'solns_spt_cvxr_drts_6*')
res_supp_6   = extract_cvxr_files(cvxr_supp_files_6,'cvxr_supp_1hr') 

cvxr_supp_files_18 = list.files(sol_dir,pattern = 'solns_spt_cvxr_drts_18*')
res_supp_18   = extract_cvxr_files(cvxr_supp_files_18,'cvxr_supp_3hr') 

# get pin2lat solutions
pin2lat_files = list.files(sol_dir,pattern='solns_pin2lat')
res_p2l = c(1:length(pin2lat_files)) %>% mclapply(function(ii){
  fname = pin2lat_files[ii]
  dat       = readRDS(paste0(sol_dir,fname))
  dat$ncp   = abs(dat$ncp)
  best_rep  = dat[which.max(dat$ncp),]$rep
  mt        = dat[dat$rep==best_rep,] %>% {.$time} %>% unlist()
  Nmeas     = length(mt)
  ncp       = costfun_svdpower(mt,freqs,Amp=NaN,cfuntype = 'ncp')
  stat_code = "NA" 
  stat_desc = NA 
  sol_name  = 'pin2lat'
  return(data.frame(Nmeas=Nmeas,ncp=ncp,stat_code=stat_code,stat_desc=stat_desc,
             solver=sol_name))
}) %>% rbindlist() %>% data.frame()

# get uniform answers
Nm_vals = c(16:48)
res_unif=NULL
for (Nm in Nm_vals){
  mt  = c(1:Nm)/Nm - 1/Nm
  ncp = costfun_svdpower(mt,freqs,Amp=NaN,cfuntype = 'ncp')
  res_unif=rbind(res_unif,data.frame(Nmeas=Nm,ncp=ncp,stat_code='NA',stat_desc=NA,
             solver='uniform'))
}

res = rbind(res_cvxr,res_supp_6,res_supp_18,res_unif,res_p2l)
df  = data.frame(res)

df %>% ggplot(aes(x=Nmeas,y=ncp,group=solver,color=solver))+
  geom_point(aes(shape=stat_code),size=2)
