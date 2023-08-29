library("tidyverse")

write_pscan_input = function(file){
  if(!(dir.exists("pscan/pscan_files"))){
    dir.create( "pscan/pscan_files" )
  }
  outfile_name = paste0(str_replace(file, ".csv", ""), "_pscan.csv")
  gene_list = read_csv(paste0(path_to_cyclops_ordering, "diff_rhythms/enrichR_files/", file))
  if (nrow(gene_list) > 1){
    gene_list = dplyr::select(gene_list, Gene_Symbols) %>%unname %>% unlist
    library(org.Hs.eg.db)
    refseq_ids <- mapIds(org.Hs.eg.db,
                         keys=gene_list,
                         column="REFSEQ",
                         keytype="SYMBOL",
                         multiVals="first")
    out = as.data.frame(refseq_ids) %>% rownames_to_column(var = "GENE_SYMBOLS")
    write.table(out, paste0( "pscan/pscan_files/", outfile_name ), sep = ',', row.names = F, col.names = T, quote = F)
  }

}

