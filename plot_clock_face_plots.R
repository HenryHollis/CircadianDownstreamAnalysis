library(NMOF)
library("JuliaCall")
julia = julia_setup("/Applications/Julia-1.6.app/Contents/Resources/julia/bin/")

cosine_distance <- function(xs, ys) {
  1 - cos(xs - ys)
}

gridfunc = function(param, true_phases, pred_phases){
  shift_from_original = param[1]
  mean(cosine_distance(true_phases, (pred_phases + shift_from_original)%%(2*pi))^2)
}

find_best_forward_backward_alignment_grid_search <- function(l1, l2) {
  true_phases = l1[!is.na(l2)]
  pred_phases = l2[!is.na(l2)]
  
  forward_search = gridSearch(gridfunc, levels = list(x =  seq(-2*pi, 2*pi, by = .05)), method = 'multicore', mc.control = list(mc.cores = 10), true_phases = true_phases, pred_phases = pred_phases)
  forward_min = forward_search$minfun
  forward_shift = forward_search$minlevel
  forward_list_of_phases = (pred_phases+forward_shift)%%(2*pi)
  
  pred_phases = (-l2[!is.na(l2)])%%(2*pi)
  
  reverse_search = gridSearch(gridfunc, levels = list(x =  seq(-2*pi, 2*pi, by = .05)), method = 'multicore', mc.control = list(mc.cores = 10), true_phases = true_phases, pred_phases = pred_phases)
  reverse_min = reverse_search$minfun
  reverse_shift = reverse_search$minlevel
  reverse_list_of_phases = (pred_phases+reverse_shift)%%(2*pi)
  
  if(forward_min < reverse_min) {
    return(forward_list_of_phases)
  }
  return(reverse_list_of_phases)
}

mouse_data = data.frame(acrophase = c(0, 0.079063, 0.151440, 2.29555, 2.9090, 2.9870,
                     2.991, 3.007, 3.12197, 3.3058, 3.31357, 3.42557,
                     3.5007, 3.8865, 4.99480, 5.04951, 6.0077),
                     Gene_Symbols = toupper(c("Arntl", "Clock", "Npas2", "Nr1d1", "Bhlhe41", "Nr1d2", 
                     "Dbp", "Ciart", "Per1", "Per3", "Tef", "Hlf", "Cry2",
                     "Per2", "Cry1", "Rorc", "Nfil3")))
mouse_data_all_cell_eneuron = mouse_data
mouse_data_all_cell_eneuron$Gene_Symbols = paste0(mouse_data_all_cell_eneuron$Gene_Symbols , "_Eneurons")
mouse_data_all_cell_astro = mouse_data
mouse_data_all_cell_astro$Gene_Symbols = paste0(mouse_data_all_cell_astro$Gene_Symbols , "_Astro")
mouse_data_all_cell = bind_rows(mouse_data_all_cell_astro, mouse_data_all_cell_eneuron)

plot_clock_face = function(plotname, df_filename,mouse_data = mouse_data, BHQ_cutoff=0.05, amp_ratio_cutoff = 0.1,best_align = F){
  df = read_csv(df_filename, show_col_types = F)
  df = dplyr::filter(df, amp_ratio >= amp_ratio_cutoff & BHQ < BHQ_cutoff)
  keep_genes = intersect(df$Gene_Symbols, mouse_data$Gene_Symbols)
  #filter out df and mouse list to just the genes in common:
  df = filter(df, Gene_Symbols %in% keep_genes) %>% arrange(Gene_Symbols)
  mouse_data = filter(mouse_data, Gene_Symbols %in% keep_genes) %>% arrange(Gene_Symbols)
  if (!best_align & any(df$Gene_Symbols == "ARNTL")){
    print("ARNTL set to 0 radians")
    subtract_acro = df$acrophase[df$Gene_Symbols == "ARNTL"]
    df$phase_MA = df$acrophase - subtract_acro
    if(sum(cosine_distance(df$phase_MA, mouse_data$acrophase)) > sum(cosine_distance(-df$phase_MA, mouse_data$acrophase))){
      df$phase_MA = -df$phase_MA
    }
  }else{
    print("ARNTL not found, using best alignment")
    df$phase_MA = find_best_forward_backward_alignment_grid_search(mouse_data$acrophase, df$acrophase)
    
  }
  
  #plot(df$phase_MA, mouse_data$acrophase)
  
  setwd("~/Box Sync/Henry_stuff/AD_project/scROSMAP/Rscripts/automatic_downstream_analysis/")
  if(length(keep_genes )>1){
    julia_source("plot_clock_face.jl")
    
    julia_call(
      "plot_clock_face",
      plotname, df$Gene_Symbols, df$phase_MA, mouse_data$Gene_Symbols, mouse_data$acrophase, df$BHQ, df$amp_ratio,
      need_return = c("R", "Julia", "None"),
      show_value = F
    ) 
  }else{
    print("Not enough genes meet cutoffs, change BHQ_cutoff or amp_ratio_cutoff")
  }
}
