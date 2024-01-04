library(tidyverse)
library(doParallel)
library(gridExtra)
library(progress)

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


plot_gene_trace = function(cyc_pred, tmm, seedlist,  useBatch = F,
                           percentile = 0.025, savePlots = F, split_cond_plots = T){
  if(useBatch){print("NOTE: Using batches in plotting")}
  
  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% arrange(Phase)
    
  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase) %>% arrange(Phase)
  }
  
   df = column_to_rownames(tmm, var = "Gene_Symbols") %>%
     t %>% as.data.frame %>% rownames_to_column( var = "Subjid") %>%
     merge(., cyc_pred, by.x = "Subjid", by.y = "ID") %>% arrange(Phase)

  if (useBatch){
    b = as.factor(df$Batch_D)      #batch factor
  }
  I = as.factor(df$Cond_D)  # condtion factor
  times = as.numeric(df$Phase) #in the case that I have CYCLOPS preds for subs not in tmm...
  
  all_genes = foreach (i = 1:length(seedlist)) %do%{
    if (!(seedlist[i] %in% colnames(df))) {print(paste(seedlist[i],"not found")); return()}
    gexp1 = as.numeric(unlist(df[,seedlist[i]]))
    times1 = df$Phase
    I1 = I
    my_df = df
    if(useBatch){b1 = b}
    
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(df))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){ #if there are NA data points...
        gexp1 = gexp1[-rm_NA] #remove them from the expression, times, and condition vectors
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        my_df = my_df[-rm_NA,]
        if(useBatch){b1 = b1[-rm_NA]} #and batch vector if using batch
      }
      
      #blunt x percentile in each condition separately
      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)
      
      if (useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + b1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1)
      }else{
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + 0)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + 0)
      }
      anova_results = anova(partial_model, full_model)
      
      if(useBatch){
        #When you have multiple batches, which batch do you use as the fitted values? We take the weighted average:
        CTL_B1 = which(b == levels(b)[1] & I1 == levels(I1)[1])
        CTL_B1_mesor = full_model[["coefficients"]][["(Intercept)"]] 
        CTL_B2 = which(b == levels(b)[2] & I1 == levels(I1)[1])
        CTL_B2_mesor = full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["b1cond_1"]] 
        CTL_either = which(I1 == levels(I1)[1])
        avg_cond_0_fitted = ( CTL_B1_mesor * length(CTL_B1) + 
                              CTL_B2_mesor * length(CTL_B2)  ) / length(CTL_either)
        full_model$fitted.values[CTL_B1] = full_model$fitted.values[CTL_B1] + (avg_cond_0_fitted - CTL_B1_mesor)
        full_model$fitted.values[CTL_B2] = full_model$fitted.values[CTL_B2] + (avg_cond_0_fitted - CTL_B2_mesor)
        
        AD_B1 = which(b == levels(b)[1] & I1 == levels(I1)[2])
        AD_B1_mesor = full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["I1cond_1"]] 
        AD_B2 = which(b == levels(b)[2] & I1 == levels(I1)[2])
        AD_B2_mesor = full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["I1cond_1"]] + full_model[["coefficients"]][["b1cond_1"]]
        AD_either = which(I1 == levels(I1)[2])
        avg_cond_1_fitted = ( AD_B1_mesor * length(AD_B1) + 
                               AD_B2_mesor * length(AD_B2) ) / length(AD_either)
        
        full_model$fitted.values[AD_B1] = full_model$fitted.values[AD_B1] + (avg_cond_1_fitted - AD_B1_mesor)
        full_model$fitted.values[AD_B2] = full_model$fitted.values[AD_B2] + (avg_cond_1_fitted - AD_B2_mesor)
      }
      
      fitted_ctl = full_model$fitted.values[I1==levels(I1)[1]]
      fitted_ad = full_model$fitted.values[I1==levels(I1)[2]]
    
      my_df[,seedlist[i]] =  gexp1
      ylim_max = max(my_df[,seedlist[i]])
      ylim_min = min(my_df[,seedlist[i]])
      df_AD = my_df %>% dplyr::filter(Cond_D == "cond_1")
      df_AD$fitted_values = fitted_ad
      df_CTL = my_df %>% dplyr::filter(Cond_D == "cond_0")
      df_CTL$fitted_values = fitted_ctl
      if(split_cond_plots){
        
        p1 = ggplot(df_CTL , aes(x = Phase , y = df_CTL[,seedlist[i]])) +
          ylim(ylim_min, ylim_max)+
          geom_point(aes(color = "CTL")) +
          geom_line(mapping=aes(x=Phase, y=fitted_values, color = "CTL"), linetype = "solid",linewidth = 2) +
          labs(title = paste0(seedlist[i], " in CTL"), x = "Circadian Phase", y = "Expression")+
          #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
          #scale_shape_manual(values=c(2, 16))+
          scale_colour_manual(values = c("blue"))
        
       
        p2 = ggplot(df_AD , aes(x = Phase , y = df_AD[,seedlist[i]])) +
          ylim(ylim_min, ylim_max)+
          geom_point(aes(color = "AD")) +
          geom_line(mapping=aes(x=Phase, y=fitted_values, color = "AD"), linetype = "solid",linewidth = 2) +
          labs(title = paste0(seedlist[i], " in AD"), x = "Circadian Phase", y = "Expression")+
          #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
          #scale_shape_manual(values=c(2, 16))+
          scale_colour_manual(values = c("red"))
        p = grid.arrange(p1, p2, nrow = 1)
      }else{
        p = ggplot(df_AD, aes(x = Phase , y = df_AD[,seedlist[i]], color = "AD")) +
          geom_point(aes(color = "AD")) +
          geom_line(mapping=aes(x=Phase, y=df_AD[,seedlist[i]], color = "AD"), linetype = "solid", linewidth = 2) +
          geom_point(data = df_CTL, mapping = aes(x = Phase , y = df_CTL[,seedlist[i]], color = "CTL")) +
          geom_line(data=df_CTL, mapping=aes(x=Phase, y=fitted_values, color = "CTL"), linetype = "solid", linewidth = 2) +
          labs(title = paste0(seedlist[i]), x = "Predicted Phase", y = "Expression")+
          scale_colour_manual(values = c("red", "blue"))
        
      }
      
      
      #print(p)
      if(savePlots){
        ggsave(paste0("plots/", seedlist[i], "_gene_trace.png"),p)
      }
    }
  }

}


plot_core_clock_genes = function(tmm_path, cyclops_path, seedlist = NULL, useBatch = useBatch,
                                 percentile = percentile, split_cond_plots = T){
  tmm = read_csv(tmm_path, show_col_types = FALSE)
  cyc_pred_file = list.files(path = paste0(cyclops_path, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(cyclops_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)
  setwd(paste0(cyclops_path, "/downstream_output"))
  if (!(dir.exists("plots"))){
    dir.create("plots")
  }
  
  genelist = c("ARNTL", "NPAS2", "CLOCK", "CRY1", "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "TEF")
  if(is.null(seedlist)){
    seedlist = genelist
  }
  plot_gene_trace(cyc_pred, tmm, seedlist, useBatch = useBatch,
                  savePlots = T, percentile = percentile, split_cond_plots = T)
  
  
}

plot_subject_histogram = function(cyclops_path, cond_subset){
  
  cyc_pred_file = list.files(path = paste0(cyclops_path, "/Fits/"), pattern = '*Fit_Output_*')
  fits = read_csv(paste(cyclops_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)
  setwd(cyclops_path)
  
  fits = filter(fits, Covariate_D == cond_subset)
  bin = 2*pi/8
  title = if(cond_subset == "cond_0") "CTL subjects" else "AD subjects"
  jpeg(file=paste0("downstream_output/plots/", cond_subset, "_phase_histogram.jpg"))
  hist(fits$Phase, breaks = seq(0, 2*pi, by = pi/6), xaxt = "n", main = title)
  axis(side=1, at=c(0,pi, 2*pi),
       labels=c("0",expression(pi),expression(2*pi)))
  dev.off()
}


