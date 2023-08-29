library(tidyverse)


write_rnks = function(path_to_cyclops_ordering){
  ## is cycling in CTL
  CTL_cyclers = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/cosinor_results_CTL.csv"))
  ranked_CTL_cyclers = arrange(CTL_cyclers, p_statistic)
  df1 = data.frame(genes = ranked_CTL_cyclers$Gene_Symbols, metric = -log(ranked_CTL_cyclers$p_statistic))
  write.table(df1, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/CTL_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  ## is cycling in AD
  AD_cyclers = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/cosinor_results_AD.csv"))
  ranked_AD_cyclers = arrange(AD_cyclers, p_statistic)
  df2 = data.frame(genes = ranked_AD_cyclers$Gene_Symbols, metric = -log(ranked_AD_cyclers$p_statistic))
  write.table(df2, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/AD_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### DR genes #####
  # AR .1, AD/CTL ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio1.csv"))
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, Log_AD_CTL_ampRatio)
  df3 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = ranked_DR_genes_AR1$Log_AD_CTL_ampRatio)
  write.table(df3, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/DRgenesAmpRatio1_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .1, -log ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio1.csv"))
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, p)
  df4 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = -log(ranked_DR_genes_AR1$p))
  write.table(df4, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/DRgenesAmpRatio1_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25, -log(p) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio25.csv"))
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, Log_AD_CTL_ampRatio)
  df5 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = ranked_DR_genes_AR25$Log_AD_CTL_ampRatio)
  write.table(df5, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/DRgenesAmpRatio25_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #AR .25, -log(p) ranked
  ranked_DR_genes_AR25 = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio25.csv"))
  ranked_DR_genes_AR25 = arrange(ranked_DR_genes_AR25, p)
  df6 = data.frame(genes = ranked_DR_genes_AR25$Gene_Symbols, metric = -log(ranked_DR_genes_AR25$p))
  write.table(df6, paste0(path_to_cyclops_ordering, "diff_rhythms/fGSEA/rnk_files/DRgenesAmpRatio25_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
}
