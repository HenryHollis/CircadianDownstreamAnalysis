library(tidyverse)
library(doParallel)
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



#test which genes are cycling from cyclops subject phase prediction
is_cycling = function(cyc_pred, tmm, cond_subset, pb = NULL, percentile = 0.025){
  print(paste("Running is_cycling() on cond_subset:", cond_subset))
  #test significant in the following genes, here that all of them.
  seedlist = unlist(unname(tmm[!grepl("_D", unlist(tmm[,1])), 1])) #ASSUMES FIRST COL is names
  cond_row_of_tmm = which(tmm[, 1] == "Cond_D")
  batch_row_of_tmm = which(tmm[, 1] == "Batch_D")
  
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% filter(Covariate_D == cond_subset) %>% arrange(Phase)
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # since seedlist is all genes, "gene" will be tmm without gene_names
  gene = apply(gene, 2, as.numeric)
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  # get the transpose, subjects x genes and put in order of CYCLOPS order
  colnames(gene1) =  unname(unlist(tmm[which(unname(unlist(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1
  
  
  b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) #the batch covariate
  #times = as.numeric(preds$Phase)
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #use "match" in the case I have CYCLOPS preds for subs not in tmm...
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    b1 = b
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        b1 = b1[-rm_NA]
      }
      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      partial_model1 = lm(gexp1 ~ b1)
      #partial_model1 = lm(gexp1 ~ 0)
      full_model1 = lm(gexp1 ~ sin(times1) + cos(times1)+ b1)
      #full_model1 = lm(gexp1 ~ sin(times) + cos(times) + b)
      anova_results1 = anova(partial_model1, full_model1)
      sin_coff = full_model1[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model1[["coefficients"]][["cos(times1)"]]
      acrophase = atan2(sin_coff, cos_coeff) %% (2*pi)
      
      p_statistic = anova_results1$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      amplitude = sqrt(sin_coff^2 + cos_coeff^2)
      #amp_ratio = amplitude/ full_model1[["coefficients"]][1]
      #When you have multiple batches, which batch do you use as the mesor. We take the weighted average
      amp_ratio = amplitude / (( full_model1[["coefficients"]][1] * length(which(b1 == levels(b1)[1])) +
                                           (full_model1[["coefficients"]][4] + full_model1[["coefficients"]][1])* length(which(b1 == levels(b1)[2])) )
                                        /length(b) ) #full_model1[["coefficients"]][1]
      #abs_amp_ratio = abs(amp_ratio)
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }
      gene_summary = cbind( Gene_Symbols, acrophase,amplitude, p_statistic, amp_ratio, sin_coff, cos_coeff)
      
      return(gene_summary)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA))
    
  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_statistic), "bonferroni")
  return(all_genes)
  
}

diff_rhyth = function(cyc_pred, tmm, seedlist,  pb = NULL, percentile = 0.025){
  print(paste("Running diff_rhyth() on seedlist of size:", length(seedlist)))
  cond_row_of_tmm = which(tmm[, 1] == "Cond_D")
  batch_row_of_tmm = which(tmm[, 1] == "Batch_D")
  
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  preds = cyc_pred %>% dplyr::select(ID, Phase, Covariate_D, batch) %>% arrange(Phase)
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1
  
  b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) 
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])        # CTL or AD factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I1 = I
    b1 = b
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        b1 = b[-rm_NA]
      }
      
      #blunt groups separately
      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)
      
      partial_model1 = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + b1)
      
      full_model1 = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1)
      
      anova_results1 = anova(partial_model1, full_model1)
      
      p = anova_results1$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      sin_coeff = full_model1[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model1[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model1[["coefficients"]][["I1cond_1:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model1[["coefficients"]][["I1cond_1:cos(times1)"]] + cos_coeff
      acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
      acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
      amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }
      
      info = cbind( Gene_Symbols, p, acrophase_AD, acrophase_CTL, amplitude_AD, amplitude_CTL)
      return(info)
    }  
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA))
  }
  
  
  all_genes = as_tibble(all_genes)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p), "bonferroni")
  all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)
  
}

diff_rhyth_AD_severity = function(cyc_pred, tmm, seedlist, rosmap_clin_path,  pb = NULL, percentile = 0.025){
  #right now this selects just AD subjects and tests the following LM:
  #y = sin(t) + cos(t) + disease_sev + batch
  #y = disease_sev*sin(t) + disease_sev*cos(t) + disease_sev + batch
  print("Running diff_rhyth_AD_severity()")
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
  cond_row_of_tmm = which(tmm[, 1] == "Cond_D")
  batch_row_of_tmm = which(tmm[, 1] == "Batch_D")
  
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  cyc_pred_merged = merge(cyc_pred, rosmap_clin, by.x = "ID", by.y = "projid", y.keep = F)
  
  preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, batch) %>% arrange(Phase)
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1
  
  cog = as.factor(preds$cogdx[match(rownames(gene1), preds$ID)])      # cogdx score 4 or 5
  cerad = as.factor(preds$ceradsc_bin[match(rownames(gene1), preds$ID)])
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) 
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I_local_cog = cog
    I_local_cerad = cerad
    b1 = b
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I_local_cog = cog[-rm_NA]
        I_local_cerad = cerad[-rm_NA]
        b1 = b1[-rm_NA]
      }
      
      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      partial_model1 = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + b1)
      full_model1 = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1)
      anova_results1 = anova(partial_model1, full_model1)
      p_cog = anova_results1$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      
      sin_coeff = full_model1[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model1[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model1[["coefficients"]][["I_local_cog5:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model1[["coefficients"]][["I_local_cog5:cos(times1)"]] + cos_coeff
      acrophase_cog4 = atan2(sin_coeff, cos_coeff) %% (2*pi)
      amplitude_cog4= sqrt((sin_coeff^2) + (cos_coeff^2))
      acrophase_cog5 = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_cog5= sqrt((sin_coeff2^2) + (cos_coeff2^2))
      
      ####### ceradsc_binned #########
      partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + b1)
      full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + b1)
      anova_results_cerad = anova(partial_model_cerad, full_model_cerad)
      p_cerad = anova_results_cerad$`Pr(>F)`[2]
      
      sin_coeff_cerad = full_model_cerad[["coefficients"]][["sin(times1)"]]
      cos_coeff_cerad = full_model_cerad[["coefficients"]][["cos(times1)"]]
      sin_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):sin(times1)"]] + sin_coeff_cerad
      cos_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):cos(times1)"]] + cos_coeff_cerad
      acrophase_cerad1to2 = atan2(sin_coeff_cerad, cos_coeff_cerad) %% (2*pi)
      amplitude_cerad1to2 = sqrt((sin_coeff_cerad^2) + (cos_coeff_cerad^2))
      acrophase_cerad3to5 = atan2(sin_coeff2_cerad, cos_coeff2_cerad) %% (2*pi)
      amplitude_cerad3to5= sqrt((sin_coeff2_cerad^2) + (cos_coeff2_cerad^2))
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }
      
      info = c( Gene_Symbols, p_cog, p_cerad, acrophase_cog4, acrophase_cog5, 
                acrophase_cerad1to2, acrophase_cerad3to5, amplitude_cog4, 
                amplitude_cog5, amplitude_cerad1to2, amplitude_cerad3to5)
      return(info)
    }  
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA, NA, NA, NA, NA))
  }
  
  
  all_genes = as_tibble(all_genes)
  colnames(all_genes) = c("Gene_Symbols", "p_cogdx","p_ceradsc", "acrophase_cog4", "acrophase_cog5", 
                          "acrophase_cerad1to2", "acrophase_cerad3to5", "amplitude_cog4", 
                          "amplitude_cog5", "amplitude_cerad1to2", "amplitude_cerad3to5")
  all_genes$BHQ_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "BH")
  all_genes$Bonf_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "bonferroni")
  all_genes$BHQ_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "BH")
  all_genes$Bonf_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "bonferroni")
  #all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)
  
}

mesor_differences = function(cyc_pred, tmm, DR_genes, pb = NULL, percentile = 0.025){ ##
  print("Running Mesor_differences()")
  cond_row_of_tmm = which(tmm[, 1] == "Cond_D")
  batch_row_of_tmm = which(tmm[, 1] == "Batch_D")
  
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  preds = cyc_pred %>% dplyr::select(ID, Phase, Covariate_D, batch) %>% arrange(Phase)
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% DR_genes), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% DR_genes), 1]))  #add the gene names to the columns of gene1
  
  b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) # sequencing batch
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])        # CTL or AD factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
    
    partial_model1 = lm(gexp1 ~ sin(times) + cos(times) + b)

    full_model1 = lm(gexp1 ~ sin(times) + cos(times) + I + b)
    anova_results1 = anova(partial_model1, full_model1)
    wilcox_test = wilcox.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    p_wilcox = wilcox_test$p.value
    
    t_test = t.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    p_ttest = t_test$p.value
    
    p_mesor = anova_results1$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    #When you have multiple batches, which batch do you use as the mesor. We take the weighted average
    # here avg_cond_0_mesor is (num_cond0_b0 * intercept + num_cond0_b1 * (intercept + b1_offset)  ) / num_Cond_0
    avg_cond_0_mesor = (full_model1[["coefficients"]][["(Intercept)"]] * length(which(b == levels(b)[1] & I == levels(I)[1])) +
        (full_model1[["coefficients"]][["(Intercept)"]] + full_model1[["coefficients"]][["bcond_[4,10)"]]) * length(which(b == levels(b)[2] & I == levels(I)[1])) ) /
        length(which(I == levels(I)[1])) 
    # here avg_cond_1_mesor is (num_cond1_b0 * (intercept + cond1_offset) + num_cond1_b1 * (intercept + cond1_offset + b1_offset)  ) / num_Cond_1
    avg_cond_1_mesor = ((full_model1[["coefficients"]][["(Intercept)"]] + full_model1[["coefficients"]][["Icond_1"]] ) * length(which(b == levels(b)[1] & I == levels(I)[2])) +
                          (full_model1[["coefficients"]][["(Intercept)"]] + full_model1[["coefficients"]][["bcond_[4,10)"]] + full_model1[["coefficients"]][["Icond_1"]] ) * length(which(b == levels(b)[2] & I == levels(I)[2])) ) /
                          length(which(I == levels(I)[2])) 

    mesor_prcnt_change = avg_cond_1_mesor / avg_cond_0_mesor
    # sin_coeff = full_model1[["coefficients"]][["sin(times)"]]
    # cos_coeff = full_model1[["coefficients"]][["cos(times)"]]
    # sin_coeff2 = full_model1[["coefficients"]][["Icond_1:sin(times)"]] + sin_coeff
    # cos_coeff2 = full_model1[["coefficients"]][["Icond_1:cos(times)"]] + cos_coeff
    # acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
    # acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
    # amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
    # amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))
    if (!is.null(pb)){
      if(!pb$finished){
        pb$tick()
      }
    }
    
    info = cbind( Gene_Symbols, p_mesor, p_wilcox, p_ttest, mesor_prcnt_change)
    return(info)
  }
  all_genes = as_tibble(all_genes)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_mesor), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_mesor), "bonferroni")
  all_genes$BHQ_wilcox = p.adjust(as.numeric(all_genes$p_wilcox), "BH")
  all_genes$BHQ_ttest = p.adjust(as.numeric(all_genes$p_ttest), "BH")
  
  return(all_genes)
}


##### main function #####

run_cycling_and_dr_analysis = function(order_path, tmm_path, isCyclingBonfCutoff = 0.05, percentile = 0.025){
  tmm = read_csv(tmm_path)      #read expression data, unordered
  colnames(tmm)[1] = "gene_names" #set first column name bc sometimes they are different
  #grab these for comparison to what cyclops finds:
  cyc_pred_file = list.files(path = paste0(order_path, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(order_path, "Fits", cyc_pred_file[1], sep = '/'))
  
  #perform nested regression on CTL data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_CTL = is_cycling(cyc_pred, tmm, cond_subset = "cond_0", pb = pb, percentile = percentile)
  #record strong cyclers in CTL, for several different amplitude cutoffs
  strong_cyclers_CTL_AR25 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.25 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  strong_cyclers_CTL_AR33 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.33 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  strong_cyclers_CTL_AR1 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.1 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  
  #perform nested regression on AD data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_AD = is_cycling(cyc_pred, tmm, cond_subset = "cond_1", pb = pb, percentile = percentile)
  #record strong cyclers in AD, for several different amplitude cutoffs
  strong_cyclers_AD_AR25 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.25 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  strong_cyclers_AD_AR33 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.33 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  strong_cyclers_AD_AR1 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.1 & as.numeric(Bonf) < isCyclingBonfCutoff) %>% arrange(as.numeric(Bonf))
  
  # We only test for diff rhythmicity if a gene cycles in AD OR CTL. Here I create those unions
  seedlist_AR25 = union(strong_cyclers_AD_AR25$Gene_Symbols, strong_cyclers_CTL_AR25$Gene_Symbols)
  seedlist_AR33 = union(strong_cyclers_AD_AR33$Gene_Symbols, strong_cyclers_CTL_AR33$Gene_Symbols)
  seedlist_AR1 = union(strong_cyclers_AD_AR1$Gene_Symbols, strong_cyclers_CTL_AR1$Gene_Symbols)
  
  ##### mesor differences ######
  gene_list_mesor =  unlist(unname(tmm[!grepl("_D", unlist(tmm[,1])), 1])) # TEST ALL genes for Mesor diff (not just cyclers)
  pb <- progress_bar$new(total = length(gene_list_mesor))
  differential_mesor = mesor_differences(cyc_pred, tmm, gene_list_mesor, pb = pb, percentile = percentile)
  Ensembl = Ensembl_dict$ENSEMBL[match(differential_mesor$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  differential_mesor = cbind(Ensembl, differential_mesor) 
  
  ####### differential rhtyhms with continuous cerad covs####
  diff_rhythms_AD_severity = diff_rhyth_AD_severity(cyc_pred, tmm,
                                                    unname(unlist(strong_cyclers_AD_AR25$Gene_Symbols)),
                                                    rosmap_clin_path = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/Meta_data/cleaned_rosmap_meta_cogdxConds.csv", 
                                                    percentile = percentile)
  
  ####### differential rhythms #####
  pb <- progress_bar$new(total = length(seedlist_AR25))
  diff_rhythms25 = diff_rhyth(cyc_pred, tmm, seedlist_AR25, pb = pb, percentile = percentile)
  Ensembl = Ensembl_dict$ENSEMBL[match(diff_rhythms25$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  diff_rhythms25 = cbind(Ensembl, diff_rhythms25) 
  
  pb <- progress_bar$new(total = length(seedlist_AR33))
  diff_rhythms33 = diff_rhyth(cyc_pred, tmm, seedlist_AR33, pb = pb, percentile = percentile)
  Ensembl = Ensembl_dict$ENSEMBL[match(diff_rhythms33$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  diff_rhythms33 = cbind(Ensembl, diff_rhythms33) 
  
  pb <- progress_bar$new(total = length(seedlist_AR1))
  diff_rhythms1 = diff_rhyth(cyc_pred, tmm, seedlist_AR1, pb = pb, percentile = percentile)
  Ensembl = Ensembl_dict$ENSEMBL[match(diff_rhythms1$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  diff_rhythms1 = cbind(Ensembl, diff_rhythms1) 
  
  #adding ENSEMBL ID to results
  Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_AD$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  cycling_in_AD = cbind(Ensembl, cycling_in_AD)
  
  Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_CTL$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  cycling_in_CTL = cbind(Ensembl, cycling_in_CTL)
  
  #Create list of strong cyclers (AR 0.25 or 0.33 and Bonf < Bonfcutoff) in CTL subjects
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cyclers_CTL_AR25$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  CTL_cyclers_AR25BonfCutoff = cbind(Ensembl, strong_cyclers_CTL_AR25) %>% dplyr::select(Ensembl , Gene_Symbols )
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cyclers_CTL_AR33$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  CTL_cyclers_AR33BonfCutoff = cbind(Ensembl, strong_cyclers_CTL_AR33) %>% dplyr::select(Ensembl , Gene_Symbols )
  
  #Create list of strong cyclers (AR 0.25 or 0.33 and Bonf < Bonfcutoff) in AD subjects
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cyclers_AD_AR25$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  AD_cyclers_AR25BonfCutoff = cbind(Ensembl, strong_cyclers_AD_AR25) %>% dplyr::select(Ensembl , Gene_Symbols )
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cyclers_AD_AR33$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  AD_cyclers_AR33BonfCutoff = cbind(Ensembl, strong_cyclers_AD_AR33) %>% dplyr::select(Ensembl , Gene_Symbols )
  
  # All genes expressed in CTL and AD:
  EnrichR_background = cycling_in_CTL %>% dplyr::select(Ensembl, Gene_Symbols)
  
  #Create list of strong DR genes (AR 0.33 or 0.25 or 0.1, and BHQ< 0.2 or BHQ< 0.1, respectively) 
  DR_cyclers_AR33BHQ2 = dplyr::filter(diff_rhythms33,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_AR25BHQ2 = dplyr::filter(diff_rhythms25,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_AR1BHQ1 = dplyr::filter(diff_rhythms1,  as.numeric(BHQ) < 0.1) %>% arrange(as.numeric(BHQ))
  
  #Create list of lost cycling DR genes:
  DR_lostAmpAD_AR33BHQ2 = filter(diff_rhythms33, BHQ < 0.2, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_lostAmpAD_AR25BHQ2 = filter(diff_rhythms25, BHQ < 0.2, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_lostAmpAD_AR1BHQ1 = filter(diff_rhythms1, BHQ < 0.1, Log_AD_CTL_ampRatio < 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  
  #Create lists of gained cycling DR genes
  DR_gainAmpAD_AR33BHQ2= filter(diff_rhythms33, BHQ < 0.2, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_gainAmpAD_AR25BHQ2= filter(diff_rhythms25, BHQ < 0.2, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_gainAmpAD_AR1BHQ1 = filter(diff_rhythms1, BHQ < 0.1, Log_AD_CTL_ampRatio > 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  
  
  if (!(dir.exists(paste(order_path, "downstream_output", sep = '/')))){
    dir.create(paste(order_path, "downstream_output", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_results", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_files", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "PSEA_files", sep = '/'))
    
  }
  #Create string for isCyclingBonfCutoff.  E.g. 0.05 -> "05"
  isCyclingBonfCutoff_str = str_extract(as.character(isCyclingBonfCutoff), "(?<=\\.)\\d+")
  
  #write out all results of cycling and DR analysis
  write.table(diff_rhythms25, paste(order_path, "downstream_output","diff_rhythms_AmpRatio25.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms33, paste(order_path, "downstream_output","diff_rhythms_AmpRatio33.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms1, paste(order_path, "downstream_output","diff_rhythms_AmpRatio1.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_CTL, paste(order_path, "downstream_output","cosinor_results_CTL.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_AD, paste(order_path, "downstream_output","cosinor_results_AD.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #gene lists for enrichR
  write.table(CTL_cyclers_AR25BonfCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR25Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(CTL_cyclers_AR33BonfCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR33Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR25BonfCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR25Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR33BonfCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR33Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  #background genes for enrichR
  write.table(EnrichR_background, paste(order_path, "downstream_output", "enrichR_files","EnrichR_background.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #DR genes for enrichR
  write.table(DR_cyclers_AR33BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_cyclers_AR33BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_AR25BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_cyclers_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_AR1BHQ1, paste(order_path, "downstream_output", "enrichR_files","DR_cyclers_AR1BHQ1.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_lostAmpAD_AR33BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_lostAmpAD_AR33BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_lostAmpAD_AR25BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_lostAmpAD_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_lostAmpAD_AR1BHQ1, paste(order_path, "downstream_output", "enrichR_files","DR_lostAmpAD_AR1BHQ1.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_gainAmpAD_AR33BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_gainAmpAD_AR33BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_gainAmpAD_AR25BHQ2, paste(order_path, "downstream_output", "enrichR_files","DR_gainAmpAD_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(DR_gainAmpAD_AR1BHQ1, paste(order_path, "downstream_output", "enrichR_files","DR_gainAmpAD_AR1BHQ1.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #Mesor differences
  write.table(differential_mesor, paste(order_path, "downstream_output", "differential_mesor_all_genes.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  sig_diff_mesor = filter(differential_mesor, as.numeric(Bonf) < 0.05 ) %>% select( Ensembl, Gene_Symbols)
  write.table(sig_diff_mesor, paste(order_path, "downstream_output","enrichR_files","diff_mesor_all_genes_Bonf05.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  
  #Continuous AD differences
  write.table(diff_rhythms_AD_severity, paste(order_path, "downstream_output", "diff_rhythms_AD_severity_AR25.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  strong_cogdx_diffs = filter(diff_rhythms_AD_severity, BHQ_cogdx< 0.1) %>% dplyr::select(Gene_Symbols)
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cogdx_diffs$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  strong_cogdx_diffs = cbind(Ensembl, strong_cogdx_diffs) %>% dplyr::select(Ensembl , Gene_Symbols )
  write.table(strong_cogdx_diffs, paste(order_path, "downstream_output", "enrichR_files", "strong_cogdx_diffs_AR25.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  
  #create lists of genes for PSEA
  PSEA_CTL_cyclers_AR25BonfCutoff = strong_cyclers_CTL_AR25 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_CTL_cyclers_AR25BonfCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_CTL_cyclers_AR25Bonf", isCyclingBonfCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  PSEA_AD_cyclers_AR25BonfCutoff = strong_cyclers_AD_AR25 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_AD_cyclers_AR25BonfCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_AD_cyclers_AR25Bonf", isCyclingBonfCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  PSEA_DR_AR25BHQ2_acrodiffs = DR_cyclers_AR25BHQ2 %>% mutate(acro_diff = (as.numeric(acrophase_AD) - as.numeric(acrophase_CTL))*12/pi ) %>%
    dplyr::select(Gene_Symbols, acro_diff)
  write.table(PSEA_DR_AR25BHQ2_acrodiffs, paste0(order_path, "/downstream_output/PSEA_files/PSEA_DR_AR25BHQ2_acrodiffs.txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  
  #write out nice summary of cycling and DR genes
  summary = data.frame(List = c("Using batch in Regression", "isCyclingBonfCutoff", "Blunting_Percentile", paste0("CTL_cyclers_AR25Bonf", isCyclingBonfCutoff_str), paste0("CTL_cyclers_AR33Bonf", isCyclingBonfCutoff_str),
      paste0("AD_cyclers_AR25Bonf", isCyclingBonfCutoff_str), paste0("AD_cyclers_AR33Bonf", isCyclingBonfCutoff_str),
      "DR_cyclers_AR1BHQ1", "DR_cyclers_AR25BHQ2", "DR_cyclers_AR33BHQ2",
      "DR_lostAmpAD_AR1BHQ1", "DR_lostAmpAD_AR25BHQ2", "DR_lostAmpAD_AR33BHQ2",
      "DR_gainAmpAD_AR1BHQ1", "DR_gainAmpAD_AR25BHQ2", "DR_gainAmpAD_AR33BHQ2"), 
      Num_genes = c(1, isCyclingBonfCutoff, percentile, dim(CTL_cyclers_AR25BonfCutoff)[1], dim(CTL_cyclers_AR33BonfCutoff)[1],
      dim(AD_cyclers_AR25BonfCutoff)[1],dim(AD_cyclers_AR33BonfCutoff)[1],
      dim(DR_cyclers_AR1BHQ1)[1],dim(DR_cyclers_AR25BHQ2)[1], dim(DR_cyclers_AR33BHQ2)[1],
      dim(DR_lostAmpAD_AR1BHQ1)[1], dim(DR_lostAmpAD_AR25BHQ2)[1],dim(DR_lostAmpAD_AR33BHQ2)[1], 
      dim(DR_gainAmpAD_AR1BHQ1)[1], dim(DR_gainAmpAD_AR25BHQ2)[1], dim(DR_gainAmpAD_AR33BHQ2)[1])
  )
  write.table(summary, paste(order_path, "downstream_output","cosinor_DR_summary.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  
  
}
