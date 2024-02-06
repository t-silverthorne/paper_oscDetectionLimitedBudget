benchmark_fdr=function(mt,Nmc,fmin,fmax,Ampvals,true_freq_vals,
                       Nfreq_regr_vals,pmethod){
  Nfreq_regr_vals %>% lapply(
    function(Nfreq_regr){
    true_freq_vals %>% lapply(
      function(true_freq){
        Ampvals %>% lapply(
          function(Amp){
            Xdat       = make_simulated_data(mt,Nmc,Amp,Amp,true_freq,true_freq)
            L          = freqsweep_regr(mt,Xdat,fmin,fmax,Nfreq_regr,return_type)
            qvals      = matrix_1d_padjust(L$pvalues,1,pmethod)
            dom_freq_q = extract_dominant_freq(qvals,L$amps,L$acros)
            pdetect_q  = (dom_freq_q$min_pq_osc <.05) %>% mean()
            
            dom_freq_p = extract_dominant_freq(L$pvalues,L$amps,L$acros)
            pdetect_p  = (dom_freq_p$min_pq_osc <.05) %>% mean()
            return(data.frame(pdetect_q=pdetect_q,
                              pdetect_p=pdetect_p,
                              Amp=Amp,
                              true_freq=true_freq,
                              Nfreq_regr=Nfreq_regr)) 
          }
        ) %>% data.table::rbindlist()
      }
    ) %>% data.table::rbindlist()
      
    }
  ) %>% data.table::rbindlist()
}


