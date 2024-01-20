benchmark_fdr=function(mt,Nmc,fmin,fmax,Ampvals,Nfreq_regr_vals){
  Nfreq_regr_vals %>% lapply(
    function(Nfreq_regr){
    true_freq_vals %>% lapply(
      function(true_freq){
        Ampvals %>% lapply(
          function(Amp){
            Xdat     = make_simulated_data(mt,Nmc,Amp,Amp,true_freq,true_freq)
            L        = freqsweep_regr(mt,Xdat,fmin,fmax,Nfreq_regr,return_type)
            qvals    = matrix_1d_padjust(L$pvalues,1,'fdr')
            dom_freq = extract_dominant_freq(qvals,L$amps,L$acros)
            pdetect  = (dom_freq$min_pq_osc <.05) %>% mean()
            return(data.frame(pdetect=pdetect,
                              Amp=Amp,
                              true_freq=true_freq,
                              Nfreq_regr=Nfreq_regr)) 
          }
        ) %>% rbindlist()
      }
    ) %>% rbindlist()
      
    }
  ) %>% rbindlist()
}


