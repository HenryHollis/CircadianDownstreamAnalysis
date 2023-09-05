

write_kegg_map_files = function(path_to_gene_file, column_out, trans_dict){
  outfile_name = basename(path_to_gene_file)
  df = read.csv(path_to_gene_file)
  df$entrez = trans_dict$entrezgene_id[ match(df$Ensembl, trans_dict$ENSEMBL)]
  
  out = cbind(df$entrez, df[,column_out])
  write.table(out, paste0("diff_rhythms/KEGG_map_diagrams/", "Entrez_", column_out, "_", outfile_name), row.names = F,col.names = F,  sep = ",")
  
}
