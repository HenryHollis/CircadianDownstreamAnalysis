
library(fgsea)
library(rstudioapi)
library(data.table)
library(stringr)
library(purrr)
run_fgsea = function(files, gene_dict, pathways, gsea_param = 1,max_size = 500, dir = "fGSEA_results"){
  sapply(files, function(file){
  file_name = paste0("../", dir, "/", str_replace(file, ".rnk", "") , ".csv")
  print(file_name)
  score_type =  ifelse(grepl("minusLogPRanked", file_name, ignore.case = T) , "pos", "std")
  print(paste("Used", score_type, "fGSEA score type"))
  plot_name = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , ".png")
  
  ranks <- read.table(file, header=F, colClasses = c("character", "numeric"))
  #here I handle remapping my gene symbols to the "chip" file from MsigDB
  ranks_merged = merge(ranks, gene_dict, by.x = "V1", by.y = "Probe.Set.ID", all.y = F, all.x = T)
  remapped_idx = which(!(is.na(ranks_merged$Gene.Symbol))) #idx of ranks_merged found in gene_dict
  ranks_merged$V1[remapped_idx] = ranks_merged$Gene.Symbol[remapped_idx] #rename first column if found
  
  if(!is_empty(which(duplicated(ranks_merged$V1)))){
    ranks_merged = ranks_merged[-which(duplicated(ranks_merged$V1)),]   #remove duplicates if any
  }
  
  ranks_vec = ranks_merged[,2]
  names(ranks_vec) = ranks_merged[,1]
  set.seed(42)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks_vec,
                    minSize  = 5,
                    #eps = 0.0,
                    gseaParam = gsea_param,
                    maxSize  = max_size,
                    scoreType = score_type,
                    nPermSimple = 10000)
  
  fgseaRes = fgseaRes %>% arrange(pval)
  fwrite(fgseaRes, file_name)
  
  topPathwaysUp <- fgseaRes[ES > 0 & !is.na(pval)][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0 & !is.na(pval)][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, topPathwaysDown)

  p = plotGseaTable(pathways[topPathways], ranks_vec, fgseaRes, 
                gseaParam=gsea_param)
  
  ggsave(plot_name, p, bg = 'white')
  
  
  #if the rnk was ranked by -log(p), only the positive enriched pathways are important:
  if (grepl("minusLogPRanked", file_name, ignore.case = T)){
    p1 = fgsea_dotplot(fgseaRes,enrichment = "pos" , my_title = paste(tools::file_path_sans_ext(file), "pos enriched pathways"))
    plot_name2 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_dotplot",".png")
    ggsave(plot_name2, p1, bg = 'white')
    
    # if the rnk was ranked by log(ad_amp/ctl_amp), we want both pathways
  }else{
    if (!purrr::is_empty(topPathwaysUp)){
      p2 = fgsea_dotplot(fgseaRes,enrichment = "pos" ,my_title = paste(tools::file_path_sans_ext(file), "pos enriched pathways"))
      plot_name2 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_posEnrichment_dotplot",".png")
      ggsave(plot_name2, p2, bg = 'white')
    }
    if (!purrr::is_empty(topPathwaysDown)){
      p3 = fgsea_dotplot(fgseaRes,enrichment = "neg" ,my_title = paste(tools::file_path_sans_ext(file), "neg enriched pathways"))
      plot_name3 = paste0("../", dir,"/plots/", str_replace(file, ".rnk", "") , "_negEnrichment_dotplot",".png")
      ggsave(plot_name3, p3, bg = 'white')
    }
  }
  
})
}

