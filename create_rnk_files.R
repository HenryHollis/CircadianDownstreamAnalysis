#File to write out .rnk files uses in GSEA, updated Nov 2, 2023 to 
#work with my differential_rhyth_generalized.R which replaces use of compareRhythms
library(tidyverse)


write_rnks = function(path_to_cyclops_ordering, isCyclingBonfCutoff_str){
  print("Creating rnk files for fGSEA.")
  ## is cycling in CTL
  CTL_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_CTL.csv"), show_col_types = FALSE)
  ranked_CTL_cyclers = arrange(CTL_cyclers, p_statistic)
  df1 = data.frame(genes = ranked_CTL_cyclers$Gene_Symbols, metric = -log(ranked_CTL_cyclers$p_statistic))
  write.table(df1, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/CTL_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  ## is cycling in AD
  AD_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_AD.csv"), show_col_types = FALSE)
  ranked_AD_cyclers = arrange(AD_cyclers, p_statistic)
  df2 = data.frame(genes = ranked_AD_cyclers$Gene_Symbols, metric = -log(ranked_AD_cyclers$p_statistic))
  write.table(df2, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### DR genes #####
  # Cycling with AR > .1 & Bonf < cutoff Cycling -> DR testing, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio1.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, Log_AD_CTL_ampRatio)
  df3 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = ranked_DR_genes_AR1$Log_AD_CTL_ampRatio)
  write.table(df3, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .1 & Bonf < cutoff Cycling -> DR testing, -log(p) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio1.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, p)
  df4 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = -log(ranked_DR_genes_AR1$p))
  write.table(df4, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .25 & Bonf < cutoff Cycling -> DR testing, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio25.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, Log_AD_CTL_ampRatio)
  df5 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = ranked_DR_genes_AR25$Log_AD_CTL_ampRatio)
  write.table(df5, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .25 & Bonf < cutoff Cycling -> DR testing, -log(p) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio25.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, p)
  df6 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = -log(ranked_DR_genes_AR25$p))
  write.table(df6, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling via Method 2 (Bonf < cutoff & AR > 0.25) -> diff rhythms, -log(p value DR) RANKED
  ranked_DR_genes_AR25_method2 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_method2_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio25.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR25_method2 = arrange(ranked_DR_genes_AR25_method2, p)
  df7 = data.frame(genes = ranked_DR_genes_AR25_method2$Gene_Symbols, metric = -log(ranked_DR_genes_AR25_method2$p))
  write.table(df7, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_minusLogPRanked_method2.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling via Method 2 (Bonf < cutoff & AR > 0.25) -> diff rhythms, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR25_method2 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_method2_CyclingBonf",isCyclingBonfCutoff_str,"AmpRatio25.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR25_method2 = arrange(ranked_DR_genes_AR25_method2, Log_AD_CTL_ampRatio)
  df8 = data.table(genes = ranked_DR_genes_AR25_method2$Gene_Symbols, metric = ranked_DR_genes_AR25_method2$Log_AD_CTL_ampRatio)
  write.table(df8, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio25_Log(AD-CTL)ranked_method2.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### differential mesor p ranked
  ranked_DM_genes = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/differential_mesor_all_genes.csv"), show_col_types = FALSE)
  ranked_DM_genes = arrange(ranked_DM_genes, p_mesor)
  df9 = data.table(genes = ranked_DM_genes$Gene_Symbols, metric = -log(ranked_DM_genes$p_mesor))
  write.table(df9, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/differential_mesor_all_genes_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
}
