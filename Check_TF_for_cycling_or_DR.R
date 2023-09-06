#####
# Tues Sept 5 2023
# Script for taking a df containing TF_names (from pscan but also enrichR results),
# and checking if they are found in DE, DE, DM, cycling files. 


augment_tf_file = function(TF_filename, deseq_filename){
  TF_file = read.csv(TF_filename)
  DEseq_toptags = read.csv(deseq_filename)
  DR_AR1 = read.csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio1.csv"))
  cycling_CTL = read.csv(paste0(path_to_cyclops_ordering, "diff_rhythms/cosinor_results_CTL.csv"))
  cycling_AD = read.csv(paste0(path_to_cyclops_ordering, "diff_rhythms/cosinor_results_AD.csv"))
  Diff_mesor = read.csv(paste0(path_to_cyclops_ordering, "diff_rhythms/differential_mesor_AR1.csv"))
    
  TF_file$cycling_in_CTL_BHQ = cycling_CTL$BHQ[match(toupper(TF_file$TF_NAME), cycling_CTL$Gene_Symbols)]
  TF_file$cycling_in_AD_BHQ = cycling_AD$BHQ[match(toupper(TF_file$TF_NAME), cycling_AD$Gene_Symbols)]
  TF_file$DR_AR1_BHQ = DR_AR1$BHQ[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)] 
  TF_file$DR_logAmpRatio = DR_AR1$Log_AD_CTL_ampRatio[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)] 
  TF_file$diff_mesor_AR1 = Diff_mesor$BHQ[match(toupper(TF_file$TF_NAME), Diff_mesor$Gene_Symbols)]
  TF_file$DEseq_DE_BHQ = DEseq_toptags$padj[match(toupper(TF_file$TF_NAME), DEseq_toptags$X)] 
  
  write.table(TF_file, TF_filename, row.names = F, col.names = T, sep = ',')
}
