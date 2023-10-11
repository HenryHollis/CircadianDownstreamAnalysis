library(tidyverse)


write_rnks = function(path_to_cyclops_ordering){
  ## is cycling in CTL
  CTL_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_CTL.csv"))
  ranked_CTL_cyclers = arrange(CTL_cyclers, p_statistic)
  df1 = data.frame(genes = ranked_CTL_cyclers$Gene_Symbols, metric = -log(ranked_CTL_cyclers$p_statistic))
  write.table(df1, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/CTL_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  ## is cycling in AD
  AD_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_AD.csv"))
  ranked_AD_cyclers = arrange(AD_cyclers, p_statistic)
  df2 = data.frame(genes = ranked_AD_cyclers$Gene_Symbols, metric = -log(ranked_AD_cyclers$p_statistic))
  write.table(df2, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### DR genes #####
  # AR .1, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio1.csv"))
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, Log_AD_CTL_ampRatio)
  df3 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = ranked_DR_genes_AR1$Log_AD_CTL_ampRatio)
  write.table(df3, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .1, -log(p) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio1.csv"))
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, p)
  df4 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = -log(ranked_DR_genes_AR1$p))
  write.table(df4, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio25.csv"))
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, Log_AD_CTL_ampRatio)
  df5 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = ranked_DR_genes_AR25$Log_AD_CTL_ampRatio)
  write.table(df5, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25, -log(p) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio25.csv"))
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, p)
  df6 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = -log(ranked_DR_genes_AR25$p))
  write.table(df6, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25 COMPARERHYTHMS -log(p) ranked
  ranked_DR_genes_AR25_compareRhythms = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio25_COMPARERHYTHMS.csv"))
  ranked_DR_genes_AR25_compareRhythms = arrange(ranked_DR_genes_AR25_compareRhythms, p_val_DR)
  df7 = data.frame(genes = ranked_DR_genes_AR25_compareRhythms$Gene_Symbols, metric = -log(ranked_DR_genes_AR25_compareRhythms$p_val_DR))
  write.table(df7, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_minusLogPRanked_COMPARERHYTHMS.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25 COMPARERHYTHMS log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR25_compareRhythms = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AmpRatio25_COMPARERHYTHMS.csv"))
  ranked_DR_genes_AR25_compareRhythms$Log_AD_CTL_ampRatio = log(ranked_DR_genes_AR25_compareRhythms$cond_1_amp/ranked_DR_genes_AR25_compareRhythms$cond_0_amp)
  ranked_DR_genes_AR25_compareRhythms = arrange(ranked_DR_genes_AR25_compareRhythms, Log_AD_CTL_ampRatio)
  df8 = data.table(genes = ranked_DR_genes_AR25_compareRhythms$Gene_Symbols, metric = ranked_DR_genes_AR25_compareRhythms$Log_AD_CTL_ampRatio)
  write.table(df8, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_Log(AD-CTL)ranked_COMPARERHYTHMS.rnk"), sep = '\t', col.names = F, row.names = F)
  
}
