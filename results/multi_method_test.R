require(devtools)
require(dplyr)
require(CVXR)
require(gurobi)
require(annmatrix)
require(ggplot2)
require(lubridate)
require(stringr)

load_all()

# current time (to nearest second)

optim_now=function(gset){
  gc()
  # output directory
  tstamp = now() %>% toString() %>%
    strsplit('\\.') %>% {.[[1]][1]} %>% 
    str_replace(' ','__')
  dir.create(paste0('results/output/Nmeas',gset$Nmeas))
  outdir =paste0('results/output/Nmeas',gset$Nmeas,'/',
                 'Nfi_',gset$Nfine,
                 '_Nfr_',gset$Nfreq,
                 '_rL1_',gset$regL1,
                 '_msa_',gset$maxit_sa,
                 '_mbf_',gset$maxit_bfgs,
                 '_nrep_',gset$nrep,
                 '_batch_',tstamp)
  dir.create(outdir)
  
  
  extract_info_from_result=function(res,tag,Nmeas=gset$Nmeas,runtime_tot,...){
    tagsp=strsplit(tag,'_')
    if(hasArg('runtime_tot')){
      # lubridate ensures that units are always in seconds
      runtime=runtime_tot %>% lubridate::as.duration() %>% as.numeric()
    }else{
      runtime=res$runtime %>% lubridate::as.duration() %>% as.numeric()
    }
    return(annmatrix(
                 x = matrix(sort(res$mtvalue),nrow=1),
              rann = data.frame(power   = 1-res$fvalue,
                                runtime = runtime,
                                tag     = tag,
                                cts     = tagsp[[1]][1],
                                lattice = tagsp[[1]][2],
                                solver  = tagsp[[1]][3]),
              cann = data.frame(time=c(1:Nmeas))))
  }
  
  
  #multi_method_opt_osc_power=function(gset,...){
  # unpack
  freqs_global  = seq(from=1,to=24,length.out=gset$Nfreq)
  Nmeas         = gset$Nmeas
  Nfreq         = gset$Nfreq
  Nfine         = gset$Nfine        
  trace_global  = gset$trace_global  
  report_global = gset$report_global
  maxit_sa      = gset$maxit_sa
  maxit_bfgs    = gset$maxit_bfgs
  Amp_global    = gset$Amp_global
  
  tau = c(1:Nfine)/Nfine -1/Nfine
  # set up controls
  ctrl_cts_arb_bfgs = list(costfun_choice='svdpower', optim_method='L-BFGS-B',
              trace=trace_global,REPORT=report_global,maxit=maxit_bfgs)
  
  ctrl_cts_arb_sa   = list(costfun_choice='svdpower',optim_method='simul_anneal',
              tfun_choice='brownian-torus',tfun_mean=0,tfun_sd=.1,
              trace=trace_global,REPORT=report_global,maxit=maxit_sa)
  
  ctrl_cts_lat_bfgs = list(costfun_choice='svdpower_2lattice',
              optim_method='L-BFGS-B',
              trace=trace_global,REPORT=report_global,maxit=maxit_bfgs)
  
  ctrl_cts_lat_sa   = list(costfun_choice = 'svdpower_2lattice',
              optim_method = 'simul_anneal',
              tfun_choice  = 'unif-with-bdry',tscale_unif_with_bdry=1/10,
              trace=trace_global,REPORT=report_global,maxit=maxit_sa)
  
  ctrl_disc_arb_sa =list(costfun_choice='svdpower_discrete',
              optim_method='simul_anneal',
              tfun_choice='single-flip',Nfine=Nfine,Nmeas=Nmeas,
              trace=trace_global,REPORT=report_global,maxit=maxit_sa)
    
  ctrl_disc_arb_cvxr = list(costfun_choice='svdpower_discrete',
              optim_method='cvxr',
              maxit=1e9,time_limit=2,MIPGapAbs=.01,
              cvxr_verbose=T,costfun_type='Linfty',
              fmin=min(freqs_global),
              fmax=max(freqs_global),Nfreq=length(freqs_global),
              lattice_cstr='none',Nfine=Nfine,Nmeas=Nmeas,
              trace=trace_global,REPORT=report_global)
  
  ctrl_disc_lat_sa = list(costfun_choice = 'svdpower_2lattice_discrete',
               optim_method   = 'simul_anneal',
               tfun_choice    = 'unif-with-bdry-discrete',tscale= 2,
               trace=trace_global,REPORT=report_global,maxit=maxit_sa)
  
  amm=NULL
  
  for (ii in c(1:gset$nrep)){
    print(paste0('BFGS: On run ',ii,' of ',gset$nrep))
    # set up initial variables
    var0_cts_arb_bfgs = runif(gset$Nmeas) 
    var0_cts_lat_bfgs = list(x0=c(shift2=runif(1),scale1=0.5,scale2=.5,shift1=0),
                             lat1 =NULL,lat2 =NULL) 
    
    
    
    # run fixed-iteration solvers
    resu_cts_arb_bfgs = opt_osc_power(dvar0   = var0_cts_arb_bfgs, 
                                      control = ctrl_cts_arb_bfgs,
                                      freqs   = freqs_global,
                                      Amp     = Amp_global)
    
    resu_cts_lat_bfgs = wrapper_sweep_lattice(opt_osc_power,
                                    Nvals    = c(gset$Nmin_2lat:gset$Nmax_2lat),
                                    Nmeas    = Nmeas,
                                    cts_flag = T,
                                    dvar0    = var0_cts_lat_bfgs,
                                    control  = ctrl_cts_lat_bfgs,
                                    freqs    = freqs_global,
                                    Amp      = Amp_global)  
    amm=rbind(amm,extract_info_from_result(resu_cts_arb_bfgs,'cts_arb_bfgs'))
    amm=rbind(amm,extract_info_from_result(
                                 res=resu_cts_lat_bfgs$res_best,
                                 tag='cts_lat_bfgs',
                                 runtime_tot=resu_cts_lat_bfgs$time_tot))
  }
  tlim = median(amm@runtime)
  
  print(paste0('CVXR: On run ',ii,' of ',1))
  ctrl_disc_arb_cvxr$time_limit=tlim
  resu_disc_arb_cvxr = opt_osc_power(dvar0  = NULL,
               control = ctrl_disc_arb_cvxr,
               freqs   = freqs_global,
               Amp     = Amp_global,
               tau     = tau)
  amm=rbind(amm,extract_info_from_result(resu_disc_arb_cvxr,'disc_arb_cvxr'))
  rm(resu_disc_arb_cvxr)
  
  
  for (ii in c(1:gset$nrep)){
    print(paste0('Annealing: On run ',ii,' of ',1))
    var0_cts_arb_sa   = runif(gset$Nmeas) 
    var0_cts_lat_sa   = list(x0=c(shift2=1/2/Nmeas,scale1=0.5,scale2=.5,shift1=0),
                             lat1 =NULL,lat2 =NULL) 
    var0_disc_arb_sa = sample(c(1:Nfine),Nmeas) 
    var0_disc_lat_sa = list(x0=c(dx1=1,dx2=1,xshift2=Nfine/2),N1=NULL,N2=NULL)
    
    # run time-limited solvers
    resu_cts_arb_sa = opt_osc_power(dvar0   = var0_cts_arb_sa,
                                    control = ctrl_cts_arb_sa,
                                    freqs   = freqs_global,
                                    Amp     = Amp_global)
    resu_cts_lat_sa = wrapper_sweep_lattice(opt_osc_power,
                                    Nvals    = c(gset$Nmin_2lat:gset$Nmax_2lat),
                                    Nmeas    = Nmeas,
                                    cts_flag = T,
                                    dvar0    = var0_cts_lat_sa,
                                    control  = ctrl_cts_lat_sa,
                                    freqs    = freqs_global,
                                    Amp      = Amp_global) 
    resu_disc_arb_sa = opt_osc_power(dvar0  = var0_disc_arb_sa,
                                    control = ctrl_disc_arb_sa,
                                    freqs   = freqs_global,
                                    Amp     = Amp_global,
                                    tau     = tau)
    resu_disc_lat_sa = wrapper_sweep_lattice(opt_osc_power,
                                    Nvals    = c(gset$Nmin_2lat:gset$Nmax_2lat),
                                    Nmeas    = Nmeas,
                                    cts_flag = F,
                                    dvar0    = var0_disc_lat_sa,
                                    control  = ctrl_disc_lat_sa,
                                    freqs    = freqs_global,
                                    Amp      = Amp_global,
                                    tau      = tau)
  
    amm=rbind(amm,extract_info_from_result(resu_cts_arb_sa,'cts_arb_sa'))
    amm=rbind(amm,extract_info_from_result(resu_disc_arb_sa,'disc_arb_sa'))
    
    amm=rbind(amm,extract_info_from_result(resu_cts_lat_sa$res_best,'cts_lat_sa',
                                           runtime_tot =resu_cts_lat_sa$time_tot))
    amm=rbind(amm,extract_info_from_result(resu_disc_lat_sa$res_best,'disc_lat_sa',
                                           runtime_tot =resu_disc_lat_sa$time_tot))
  }
  saveRDS(amm,paste0(outdir,'amm_all.RDS'))
  saveRDS(gset,paste0(outdir,'global_settings.RDS'))
}

# global settings
gset = list(
  Nmeas             = 16,
  Nfreq             = 2^6,
  Nfine             = 288,
  trace_global      = 0,
  report_global     = 1,
  maxit_sa          = 1e3,
  maxit_bfgs        = 10, 
  Amp_global        = 2,
  timelimit_by_bfgs = T,
  Nmin_2lat         = 4,
  Nmax_2lat         = 12,
  nrep              = 100,
  regL1             = 0
)
#for (Nmeas in seq(16,24,2)){
#  gset$Nmeas=Nmeas 
#  optim_now(gset)
#}

gset$Nfreq=2^8
for (Nmeas in c(16,24)){
  gset$Nmeas=Nmeas 
  optim_now(gset)
}

#make_plot=F
#if (make_plot){
#  theme_set(theme_bw())
#  amm@'' %>% ggplot(aes(x=power,y=runtime))+
#    geom_point(aes(shape=factor(lattice,c('arb','lat'),c('none','2-lattice')),
#                   color=factor(solver,c('cvxr','bfgs','sa'),c('DCP','BFGS','annealing'))),size=2)+
#    labs(shape='Lattice constraint',
#         color='Solver',
#          y='runtime [s]',
#          x='worst case power')+
#    facet_wrap(~factor(cts,c('cts','disc'),c('Continuous time','Discretised time (5 minute intervals)')),nrow=2)+
#    xlim(c(0,1))+
#    scale_color_manual(values=c('DCP'='red','BFGS'='blue','annealing'='cyan'))
#}
