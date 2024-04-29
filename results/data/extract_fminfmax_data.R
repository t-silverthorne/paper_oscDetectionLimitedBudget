require(devtools)
require(ggplot2)
require(annmatrix)
require(ggplotify)
require(patchwork)
require(parallel)
require(data.table)
require(stringr)
require(dplyr)
require(latex2exp)
require(annmatrix)
devtools::load_all()
sol_dir   = 'results/data/fminfmax//'
sol_files = list.files(sol_dir,pattern='sol_gur')

Ampglob    = 1 
Nfreq      = 47*8 

dfreq      = .05 
Nfine      = 144
tau        = c(1:Nfine)/Nfine-1/Nfine
gur_thresh = .9

fdf = NULL
for (ii in c(1:length(sol_files))){
  fname = sol_files[ii]
  gres  = readRDS(paste0(sol_dir,fname))
 
  # get settings values from parameters 
  spart   = strsplit(fname,'.RDS') %>% {.[[1]]}
  spart   = strsplit(spart,'_')
  f_df    = spart[[1]][which(spart[[1]]=='df')+1] %>% as.numeric()  
  f_fmin  = spart[[1]][which(spart[[1]]=='fmin')+1] %>% as.numeric()  
  f_Nmeas = spart[[1]][which(spart[[1]]=='Nmeas')+1] %>% as.numeric()  

  Nfine   = length(gres$x)-1
  
  if (Nfine==144){
    mtloc      = tau[gres$x[1:Nfine]>gur_thresh]
    Nmeas      = length(mtloc)
    Nmeas_good = Nmeas==f_Nmeas
    freqs_loc  = seq(f_fmin,f_fmin+f_df,dfreq)
    pwr        = costfun_svdpower(mtloc,freqs_loc,Amp=Ampglob,cfuntype = 'power')
    
    mt_unif    = c(1:Nmeas)/Nmeas - 1/Nmeas
    pwr_unif   = costfun_svdpower(mt_unif,freqs_loc,Amp=Ampglob,cfuntype='power')
    stat_code  = gres$status
   
    
    dfloc = data.frame(
      Nmeas       = Nmeas,
      fmin        = f_fmin,
      df          = f_df,
      d_power     = pwr-pwr_unif,
      Nmeas_match = Nmeas_good
      )
    
    fdf=rbind(fdf,dfloc)
  }
}
