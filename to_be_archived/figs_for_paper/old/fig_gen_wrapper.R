files_to_run=c('results/output/Nmeas16/Nfi_288_Nfr_64_rL1_0_msa_1000_mbf_10_nrep_100_batch_2024-02-07__17:51:30amm_all.RDS',
               'results/output/Nmeas18/Nfi_288_Nfr_64_rL1_0_msa_1000_mbf_10_nrep_100_batch_2024-02-07__20:57:25amm_all.RDS')

for (f in files_to_run){
  print(paste0('Running: ',f))
  write(f,file='figs_for_paper/fnames.csv')
  rmarkdown::render('figs_for_paper/makeFig1Fig2.Rmd',knit_root_dir = getwd())
}
