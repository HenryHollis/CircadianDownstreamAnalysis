library(tidyverse)
library(doParallel)
library(gridExtra)

blunt_outliers = function(vec, percentile = 0.025){
  num =length(which(!is.na(vec)))
  blunt_n_points = round(percentile * num, 0)
  ord = sort(vec)
  upper_val = ord[num -blunt_n_points]
  lower_val = ord[blunt_n_points+1]
  
  vec[which(vec > upper_val)] = upper_val
  vec[which(vec < lower_val)] = lower_val
  return(vec)
}




draw_gene_tracings = function(cyc_pred, tmm, seedlist, savePlots = F, split_cond_plots = T, percentile = 0.025){
 
  #cyc_pred$pmi = rosmap_meta$pmi[match(cyc_pred$ID, rosmap_meta$projid)]
  #mouse_path = list.files(fits_path, pattern = "Mouse_Atlas_*")
  #mouse_aligned = read_csv(paste0(fits_path, mouse_path))
  colnames(tmm)[1] = "gene_names" #set first column name bc sometimes they are different
  cyc_pred$Covariate_D = tmm[1, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  #seedlist = unname(unlist(read.csv(seedlist)))
  #sig_cyclers = dplyr::filter(mouse_aligned, BHQ_Statistic < 0.05) %>% dplyr::select( Gene_Symbols)
  #seedlist = intersect(unlist(sig_cyclers), seedlist)
  #print(seedlist)
  
  preds_AD = dplyr::filter(cyc_pred, Covariate_D == 'cond_1') %>% dplyr::select(ID, Phase, Covariate_D) %>% arrange(Phase)
  preds_N = dplyr::filter(cyc_pred, Covariate_D == 'cond_0') %>% dplyr::select(ID, Phase, Covariate_D) %>% arrange(Phase)

  
  gene = tmm[which(tmm$gene_names %in% seedlist), -1]
  #gene = gene[,colnames(gene) %in% colnames(true_times)]
  gene1 = t(gene[,na.exclude(match(preds_AD$ID, colnames(gene)))])
  gene2 = t(gene[,na.exclude(match(preds_N$ID, colnames(gene)))])
  
  colnames(gene1)= colnames(gene2) = unname(unlist(tmm[which(tmm$gene_names %in% seedlist), 1]))
  times_AD = as.numeric(preds_AD$Phase[match(rownames(gene1), preds_AD$ID)])
  times_N = as.numeric(preds_N$Phase[match(rownames(gene2), preds_N$ID)])
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    species = colnames(gene1)[gene_i]
    print("********INFO: assuming ONE batch *********")
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times_AD
    rm_gexp1 = which(is.na(gexp1))
    if(!is_empty(rm_gexp1)){
      gexp1 = gexp1[-rm_gexp1]
      times1 = times1[-rm_gexp1]
    }
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
 
    #pmi1 = as.factor(preds_AD$pmi < 20)
    partial_model1 = lm(gexp1 ~ 1)
    full_model1 = lm(gexp1 ~ sin(times1) + cos(times1) + 1)
    anova_results1 = anova(partial_model1, full_model1)
    line_type1 = ifelse(anova_results1$`Pr(>F)`[2] < 0.05, "solid", "dashed")
    predicted_values1 = full_model1$fitted.values
    plot_df1 = data.frame(gene_expression = gexp1,  x_vals = times1, predicted_vals = predicted_values1)
    plot_df1$study1= as.factor("AD")
    
    
    gexp2 = as.numeric(unlist(gene2[,gene_i]))
    times2 = times_N
    rm_gexp2 = which(is.na(gexp2))
    if(!is_empty(rm_gexp2)){
      gexp2 = gexp2[-rm_gexp2]
      times2 = times2[-rm_gexp2]
    }
    gexp2 = blunt_outliers(gexp2, percentile = percentile)
    #pmi2 = as.factor(preds_N$pmi < 20)
    partial_model2 = lm(gexp2 ~ 1)
    full_model2 = lm(gexp2 ~ sin(times2) + cos(times2) + 1)
    anova_results2 = anova(partial_model2, full_model2)
    line_type2 = ifelse(anova_results2$`Pr(>F)`[2] < 0.05, "solid", "dashed")
    predicted_values2 = full_model2$fitted.values
   
    plot_df2 = data.frame(gene_expression = gexp2, x_vals = times2, predicted_vals = predicted_values2)
    plot_df2$study2 = as.factor("CTL")
    ylim_max = max(c(gexp1, gexp2))
    ylim_min = min(c(gexp1, gexp2))
 
    if(split_cond_plots){
      p1 = ggplot(plot_df2, aes(x = x_vals , y = gene_expression)) +
        ylim(ylim_min, ylim_max)+
        geom_point(aes(color = "CTL")) +
        geom_line(data=plot_df2, mapping=aes(x=x_vals, y=predicted_vals, color = "CTL"), linetype = line_type2,linewidth = 2) +
        labs(title = paste0(species, " in CTL"), x = "Circadian Phase", y = "Expression")+
        #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
        #scale_shape_manual(values=c(2, 16))+
        scale_colour_manual(values = c("blue"))
      
      
      p2 = ggplot(plot_df1, aes(x = x_vals , y = gene_expression))+
        ylim(ylim_min, ylim_max)+
        geom_point(aes(color = "AD")) +
        geom_line(data=plot_df1, mapping=aes(x=x_vals, y=predicted_vals, color = "AD"), linetype = line_type1, linewidth = 2) +
        labs(title = paste0(species, " in AD"), x = "Predicted Phase", y = "Expression")+
        #annotate("text", x = 3, y = 0, label = paste("DR FDR", DR_FDR))+
        #scale_shape_manual(values=c(2, 16))+
        scale_colour_manual(values = c("red"))
      p = grid.arrange(p1, p2, nrow = 1)
    }else{
      p = ggplot(plot_df1, aes(x = x_vals , y = gene_expression, color = "AD")) +
        geom_point(aes(color = "AD")) +
        geom_line(data=plot_df1, mapping=aes(x=x_vals, y=predicted_vals, color = "AD"), linetype = line_type1, linewidth = 2) +
        geom_point(data = plot_df2, mapping = aes(x = x_vals , y = gene_expression, color = "CTL")) +
        geom_line(data=plot_df2, mapping=aes(x=x_vals, y=predicted_vals, color = "CTL"), linetype = line_type2, linewidth = 2) +
        labs(title = paste0(species), x = "Predicted Phase", y = "Expression")+
        scale_colour_manual(values = c("red", "blue"))
      
    }
    
    print(p)
    if(savePlots){
       #CODE TO SAVE
      }
  }
 

}
plot_exc_neuron_genes = function(seedlist, split_cond_plots = T, percentile = 0.025,
          tmm_path = "EdgeR_filt_normed/ExcNeurons_FiltByExprDefault_TMM_combatSeqAdjusted.csv",
          fits_path = "../../../training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/TMMs_w_batch/Exc3and5_FiltByEdgeRDefault_ipBulkChenZhang_condAndBatchCovs_3EGdefault_noTransferFit/Fits/"){
  setwd("~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/")
  rosmap_meta = read.csv("../../../../../scROSMAP/Meta_data/cleaned_rosmap_meta_cogdxConds.csv")
  print(paste("Plotting from", tmm_path))
  print(paste("Using", fits_path, "Fits"))
  # fits_path = "../../../training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/Exc_Neurons_CellFiltered10Percent_ErikChenZhang_condCovs_3EG_Jun12Redo/Fits/"
  #tmm_path = "ExcNeurons_cellsFiltered1count10percentCells.csv"
  tmm = read_csv(tmm_path)
  cyc_pred_file = list.files(path = fits_path, pattern = '*Fit_Output_[0-9]+')
  cyc_pred = read.csv(paste0(fits_path,cyc_pred_file))
  
  draw_gene_tracings(cyc_pred, tmm, seedlist, split_cond_plots = split_cond_plots , percentile = percentile)
}

draw_gene_tracings_AD_continuous = function(cyc_pred, tmm, seedlist, savePlots = F,
                                            rosmap_clin_path, percentile = .025){
  
  ##### read in ROSMAP clin ####
  rosmap_clin = read_csv(rosmap_clin_path)
  rosmap_clin = rosmap_clin[ na.exclude(match(cyc_pred$ID, rosmap_clin$projid)),]
  rosmap_clin = rosmap_clin %>%
    mutate(braaksc_bin = cut(braaksc, c(0, 3, 5, 7), right = F))
  rosmap_clin = rosmap_clin %>%
    mutate(ceradsc_bin = cut(ceradsc, c(1, 3, 5), right = F))
  rosmap_clin$apoe_ordinal  = 1
  rosmap_clin$apoe_ordinal[rosmap_clin$apoe_genotype == 34 | rosmap_clin$apoe_genotype == 24] = 2
  rosmap_clin$apoe_ordinal[rosmap_clin$apoe_genotype == 44 ] = 3
  rosmap_clin$apoe_ordinal[is.na(rosmap_clin$apoe_genotype) ] = NA
  ###############
  colnames(tmm)[1] = "gene_names" #set first column name bc sometimes they are different
  
  cyc_pred_merged = merge(cyc_pred, rosmap_clin, by.x = "ID", by.y = "projid", y.keep = F)
  
  
  preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, 'cogdx') %>% arrange(Phase)
  
  gene = tmm[which(tmm$gene_names %in% seedlist), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(tmm$gene_names %in% seedlist), 1]))  #add the gene names to the columns of gene1
  
  I = as.factor(preds$cogdx[match(rownames(gene1), preds$ID)])      # CTL or AD factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    species = colnames(gene1)[gene_i]
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I_local = I
    rm_NA = which(is.na(gexp1))
    if(!is_empty(rm_NA)){
      gexp1 = gexp1[-rm_NA]
      times1 = times1[-rm_NA]
      I_local = I[-rm_NA]
    }
      
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
    partial_model1 = lm(gexp1 ~ sin(times1) + cos(times1) + I_local + 0)
    full_model1 = lm(gexp1 ~ I_local*sin(times1) + I_local*cos(times1) + I_local + 0)
    anova_results1 = anova(partial_model1, full_model1)
      
    line_type1 = ifelse(anova_results1$`Pr(>F)`[2] < 0.05, "solid", "dashed")
    predicted_values1 = full_model1$fitted.values
    plot_df1 = data.frame(gene_expression = gexp1,  x_vals = times1, predicted_vals = predicted_values1, covariate = as.factor(preds$cogdx))

    
    p = ggplot(plot_df1, aes(x = x_vals , y = gene_expression)) +
      geom_point(aes(color = covariate)) +
      geom_line(data=filter(plot_df1, covariate == levels(plot_df1$covariate)[1]), mapping=aes(x=x_vals, y=predicted_vals), color = "blue", linetype = line_type1,linewidth = 2) +
      geom_line(data=filter(plot_df1, covariate == levels(plot_df1$covariate)[2]), mapping=aes(x=x_vals, y=predicted_vals),color = "red", linetype = line_type1,linewidth = 2) +
      labs(title = paste0(species, " in AD"), x = "Circadian Phase", y = "Expression")+
      #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
      #scale_shape_manual(values=c(2, 16))+
      scale_colour_manual(values = c("blue", "red"))
      
    
  
    
    print(p)
    if(savePlots){
      #CODE TO SAVE
    }
  }
  
  
}
