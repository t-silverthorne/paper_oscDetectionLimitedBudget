helper_update_amm=function(amm,Lnow,gset,tag,wrapped=F){
  for (ii in c(1:gset$nrep)){
    res = Lnow[[ii]]
    if (wrapped){
      amm = rbind(amm,helper_extract_info_from_result(res=res$res_best,tag=tag,
                                               runtime_tot=res$time_tot)) 
    }else{
      amm = rbind(amm,helper_extract_info_from_result(res=res,tag=tag)) 
    }
  }
  return(amm)  
}