## Not using this file anymore, I write out the lost and gained files in differential_rhyth_single_cov...csv


create_DRgenelists = function(DR_path, BHQ_cutoff = 0.2){
  ### Function takes file name of file from DR cosinor regression
  ### opens it and writes out the genes that gained or lost rhythmicity 
  
  df = read_csv(DR_path)
  lost_amp = filter(df, BHQ < BHQ_cutoff, Log_AD_CTL_ampRatio < 0) %>% select(Ensembl, Gene_Symbols)
  gained_amp = filter(df, BHQ < BHQ_cutoff, Log_AD_CTL_ampRatio > 0) %>% select(Ensembl, Gene_Symbols)
  
  str1 = paste0("DR_lostAmpAD_AR", str_extract(DR_path, "\\d+"), "BHQ", str_replace((BHQ_cutoff), "\\d+\\.", ""), ".csv")
  write.table(lost_amp, str1, sep = ',', row.names = F, col.names = T)
  str2 = paste0("DR_gainAmpAD_AR", str_extract(DR_path, "\\d+"), "BHQ", str_replace((BHQ_cutoff), "\\d+\\.", ""), ".csv")
  write.table(gained_amp, str2, sep = ',', row.names = F, col.names = T)
  
}
