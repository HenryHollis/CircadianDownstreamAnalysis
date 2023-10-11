library(ggplot2)

fgsea_dotplot <- function(object, show = 15, 
                          font.size=12, my_title = "GSEA", enrichment = "all", BHQ_cutoff = 0.2) {
  object = as_tibble(object)
  if (enrichment == "pos"){
    print("Including only pos enriched pathways")
    object = dplyr::filter(object, NES > 0)
  }else if (enrichment == "neg"){
    print("Including only neg enriched pathways")
    object = dplyr::filter(object, NES < 0)
  }else{
    print("Using positive and negatively enriched pathways")
  }
  data = object %>% top_n(-1*show, wt = padj)%>% arrange(padj) %>%#-1 means "bottom" n pvalues
    dplyr::filter(padj < BHQ_cutoff)
  #fixed_labels <- str_replace_all(data$pathway, "_", " ")
  data$fixed_labels = str_replace_all(data$pathway, "_", " ")
  ggplot(data, aes(x=reorder(fixed_labels,padj, decreasing = T), y= abs(NES), color = padj, size = size)) +
    geom_point() +
    coord_flip()+
    scale_color_continuous(low="red", high="blue", name = "BHQ",
                           guide=guide_colorbar(reverse=T)) +
    #scale_size_continuous( name = "Size", range  =c(3,8)) +
    scale_x_discrete(name = "Pathway", labels = function(x) str_wrap(x, width = 30)) +
    ylab("abs(NES)") +
    labs(title = str_wrap(my_title, 50)) #  # theme_dose(12)
  
}
