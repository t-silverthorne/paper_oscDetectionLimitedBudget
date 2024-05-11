source('figs_for_paper/Fig1_aug_vs_cvxr/settings.R')
Nfine=gset$Nfine
mc_cores=Sys.getenv("SLURM_CPUS_PER_TASK") 

tau = c(1:Nfine)/Nfine -1/Nfine

mtdf = NULL
for (Nmeas in Nmeasvals){
  print(paste0('Running Nmeas=',Nmeas))
  gset$Nmeas=Nmeas
  ctrl_disc_arb_cvxr = list(costfun_choice='svdpower_discrete',
                            optim_method='cvxr',
                            maxit=1e9,time_limit=gset$timelimit_cvxr,MIPGapAbs=.01,
                            cvxr_verbose=T,costfun_type='Linfty',
                            fmin=min(freqs),
                            fmax=max(freqs),Nfreq=gset$Nfreq,
                            lattice_cstr='none',Nfine=gset$Nfine,Nmeas=gset$Nmeas,
                            trace=gset$trace_global,REPORT=gset$report_global)
  
  res=opt_osc_power(dvar0  = NULL,
                    control = ctrl_disc_arb_cvxr,
                    freqs   = freqs,
                    Amp     = NaN,
                    tau     = tau,
                    cfuntype='ncp')
  mtdf=rbind(mtdf,data.frame(time=res$mtvalue,solver='cvxr',Nmeas=Nmeas,
                             ncp=-res$fvalue))
  rm(res)
}

mtdf_aug=NULL
for (Nmeas in Nmeasvals){
  print(paste0('Running Nmeas=',Nmeas))
  gset$Nmeas=Nmeas
  gset$Nmin_2lat=floor(Nmeas/4)
  gset$Nmax_2lat=3*floor(Nmeas/4)
  dvar0=list(N1=floor(Nmeas/4),
             N2=3*floor(Nmeas/4),
             shift2=1/4/Nmeas,
             scale2=1)
  control=list(N1min=gset$Nmin_2lat,N1max=gset$Nmax_2lat,
               trace=1,REPORT=1,maxit=gset$maxit_sa,
               costfun_choice='auglattice')
  #res=opt_osc_power(dvar0=dvar0,freqs=freqs,control=control,
  #           cfuntype='ncp',
  #           regFder=10)
  res=opt_osc_power(dvar0=dvar0,freqs=freqs,control=control,
                    cfuntype='ncp')
  mtdf_aug=rbind(mtdf_aug,data.frame(time=res$mtvalue,solver='pin2lat',Nmeas=Nmeas,
                                     ncp=-res$fvalue))
}
mtdf_aug
mtdf_full=rbind(mtdf,mtdf_aug)

saveRDS(mtdf_full,'figs_for_paper/Fig1_aug_vs_cvxr/solns.RDS')