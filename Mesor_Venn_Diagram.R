library(VennDiagram)


create_mesor_venn_diag = function(path_to_cyclops_ordering, BHQcutoff = 0.05){
setwd(paste(path_to_cyclops_ordering, "diff_rhythms", sep = "/"))

mesor_file = read_csv("differential_mesor_AR1.csv")
diff_mesor_BHQ05 = filter(mesor_file, BHQ < BHQcutoff) %>% dplyr::select(Gene_Symbols) %>% unname %>% unlist
diff_meanWilcox_BHQ05 = filter(mesor_file, BHQ_wilcox < BHQcutoff) %>% dplyr::select(Gene_Symbols) %>% unname %>% unlist
diff_meanttest_BHQ05 = filter(mesor_file, BHQ_ttest < BHQcutoff) %>% dplyr::select(Gene_Symbols) %>% unname %>% unlist

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.plot =  venn.diagram(
  x = list(diff_mesor_BHQ05, diff_meanWilcox_BHQ05, diff_meanttest_BHQ05),
  category.names = c("Diff Mesor" , "Diff Mean - Wilcox" , "Diff Mean - Ttest"),
  fill = c("red", "green", "blue"),
  alpha = c(0.5, 0.5, 0.5),
  cex = 1,
  cat.fontface = 2,
  lty =1, 
  cat.just = list(c(0, 0) , c(.8,0) , c(.7,0)),
  print.mode=c("raw","percent"),
  filename = NULL,
  output=TRUE
)

jpeg(filename = paste0('mesor_differences_venn_diagramm_BHQ', str_extract(as.character(BHQcutoff), "(?<=\\.)\\d+"), ".jpg"))
grid.draw(venn.plot)
dev.off()

}
