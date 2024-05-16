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
sol_dir   = 'results/data/multi_highres/'
sol_files = list.files(sol_dir,pattern='sol_gur')

Nfreq = 47*8 
freqs = seq(1,24,length.out=Nfreq)
Nfine=144
tau = c(1:Nfine)/Nfine-1/Nfine
gur_thresh = .9

sols = NULL
for (ii in c(1:length(sol_files))){
  fname = sol_files[ii]
  gres=readRDS(paste0(sol_dir,fname))
 
  # get settings values from parameters 
  spart   = strsplit(fname,'.RDS') %>% {.[[1]]}
  spart   = strsplit(spart,'_')
  f_wreg  = spart[[1]][which(spart[[1]]=='wreg')+1] %>% as.numeric()  
  f_drts  = spart[[1]][which(spart[[1]]=='drts')+1] %>% as.numeric()  
  f_Nmeas = spart[[1]][which(spart[[1]]=='Nmeas')+1] %>% as.numeric()  

  Nfine = length(gres$x)-1
  
  if (Nfine==144){
    mtloc = tau[gres$x[1:Nfine]>gur_thresh]
  
    Nmeas = length(mtloc)
    Nmeas_good = Nmeas==f_Nmeas
    
    ncp = evalWorstNCP(mt=mtloc,freqs = freqs,Amp=1)
    acro_freq_der = evalAcroFreqGradL2(tvec=mtloc,
                                       pars=list(fmin=min(freqs),fmax=max(freqs)),
                                       Amp=1)
    stat_code = gres$status
   
    rowdat=data.frame(
      Nmeas       = Nmeas,
      wreg        = f_wreg,
      drts        = f_drts,
      ncp         = ncp,
      sc_afder    = sqrt(acro_freq_der)/Nmeas,
      stat_code   = stat_code,
      Nmeas_match = Nmeas_good
      )
    mtloc=c(mtloc,rep(Inf,48-Nmeas))
    mtloc
    am=annmatrix(matrix(mtloc,nrow=1),
                 rowdat,
                 data.frame(real_meas=mtloc<Inf))
    sols=rbind(sols,am)
  }
}

#sols_plt = sols[sols@wreg %in% c(0,1),]
#sols_plt@''%>% 
#  ggplot(aes(x=Nmeas,y=ncp,
#             shape=as.factor(wreg),color=as.factor(drts)))+geom_point(size=2)+facet_wrap(~stat_code)