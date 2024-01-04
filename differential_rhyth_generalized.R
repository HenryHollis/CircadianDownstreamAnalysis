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
is_cycling = function(cyc_pred, tmm, cond_subset, pb = NULL, useBatch = F, percentile = 0.025){
  cat(paste("\nRunning is_cycling() on cond_subset:", cond_subset))
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  #test significant in the following genes, here that all of them.
  seedlist = unlist(unname(tmm[!grepl("_D", unlist(tmm[,1])), 1])) #ASSUMES FIRST COL is names

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist

  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% filter(Covariate_D == cond_subset) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase) %>% filter(Covariate_D == cond_subset) %>% arrange(Phase)
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # since seedlist is all genes, "gene" will be tmm without gene_names
  gene = apply(gene, 2, as.numeric)
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  # get the transpose, subjects x genes and put in order of CYCLOPS order
  colnames(gene1) =  unname(unlist(tmm[which(unname(unlist(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1

  #below 2 lines I use "match" in case I am given phases for subjects not in the data
  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) #the batch variable
  }
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #the given phase of each subject

  #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }
      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      if (useBatch){
        partial_model = lm(gexp1 ~ b1)
        full_model = lm(gexp1 ~ sin(times1) + cos(times1)+ b1)
      }else{
        partial_model = lm(gexp1 ~ 1)
        full_model = lm(gexp1 ~ sin(times1) + cos(times1))
      }

      anova_results = anova(partial_model, full_model)
      sin_coff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      acrophase = atan2(sin_coff, cos_coeff) %% (2*pi)

      p_statistic = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      amplitude = sqrt(sin_coff^2 + cos_coeff^2)
      if(useBatch){
        amp_ratio = amplitude / (( full_model[["coefficients"]][["(Intercept)"]] * length(which(b1 == levels(b1)[1])) +
                                     (full_model[["coefficients"]][["b1cond_1"]] + full_model[["coefficients"]][["(Intercept)"]])* length(which(b1 == levels(b1)[2])) ) /
                                   length(b) )
      }else{
        amp_ratio = amplitude/ full_model[["coefficients"]][["(Intercept)"]]
      }

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

is_cycling_method2 = function(cyc_pred, tmm, pb = NULL, useBatch = F, percentile = 0.025){
  cat("\nRunning is_cycling_method2() (analogue of compareRhythms) on all subjects")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  seedlist = unlist(unname(tmm[!grepl("_D", unlist(tmm[,1])), 1])) #ASSUMES FIRST COL is names

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist

  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase) %>% arrange(Phase)
  }


  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # since seedlist is all genes, "gene" will be tmm without gene_names
  gene = apply(gene, 2, as.numeric)
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  # get the transpose, subjects x genes and put in order of CYCLOPS order
  colnames(gene1) =  unname(unlist(tmm[which(unname(unlist(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1

  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)])      #batch factor
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)]) # CTL or AD factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])  #in the case that I have CYCLOPS preds for subs not in tmm...

  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I1 = I
    if(useBatch){b1 = b}
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)

      if (useBatch){
        partial_model = lm(gexp1 ~ I1 + b1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1)
      }else{
        partial_model = lm(gexp1 ~ I1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 )
      }

      anova_results = anova(partial_model, full_model)

      p = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I1cond_1:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I1cond_1:cos(times1)"]] + cos_coeff
      acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
      acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
      amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))

      if(useBatch){
        #When you have multiple batches, which batch do you use as the mesor? We take the weighted average
        # here avg_cond_0_mesor is ((#_cond0_b0 * intercept) + num_cond0_b1 * (intercept + b1_offset)  ) / num_Cond_0

        avg_cond_0_mesor = ( (full_model[["coefficients"]][["(Intercept)"]] * length(which(b == levels(b)[1] & I == levels(I)[1])) ) +
                               ((full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["b1cond_1"]]) * length(which(b == levels(b)[2] & I == levels(I)[1])) ) ) /
          length(which(I == levels(I)[1]))

        # here avg_cond_1_mesor is (#_cond1_b0 * (intercept + cond1_offset) + #_cond1_b1 * (intercept + cond1_offset + b1_offset)  ) / num_Cond_1
        avg_cond_1_mesor = ((full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["I1cond_1"]] ) * length(which(b == levels(b)[1] & I == levels(I)[2])) +
                              (full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["I1cond_1"]] + full_model[["coefficients"]][["b1cond_1"]] ) * length(which(b == levels(b)[2] & I == levels(I)[2])) ) /
          length(which(I == levels(I)[2]))

        mesor_CTL = avg_cond_0_mesor
        mesor_AD = avg_cond_1_mesor

      }else{
        mesor_CTL = full_model[["coefficients"]][["(Intercept)"]]
        mesor_AD = full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["I1cond_1"]]
      }
      amp_ratio_CTL = amplitude_CTL/ mesor_CTL
      amp_ratio_AD = amplitude_AD/ mesor_AD

      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }

      info = cbind( Gene_Symbols, p, acrophase_AD, acrophase_CTL, amplitude_AD, amplitude_CTL, amp_ratio_CTL, amp_ratio_AD, mesor_CTL, mesor_AD)
      return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA, NA, NA, NA))
  }


  all_genes = as_tibble(all_genes)
  all_genes$p = as.numeric(all_genes$p)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p), "bonferroni")
  all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)


}

diff_rhyth = function(cyc_pred, tmm, seedlist,  pb = NULL, useBatch = F, percentile = 0.025){
  cat(paste("\nRunning diff_rhyth() on seedlist of size:", length(seedlist)))
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist

  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase) %>% arrange(Phase)
  }


  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1

  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)])      #batch factor
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])  # condtion factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...

  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I1 = I
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

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

      p = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I1cond_1:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I1cond_1:cos(times1)"]] + cos_coeff
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
  all_genes$p = as.numeric(all_genes$p)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p), "bonferroni")
  all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)

}

diff_rhyth_AD_severity = function(cyc_pred, tmm, seedlist, rosmap_clin_path,  pb = NULL, useBatch = F, percentile = 0.025){
  cat("\nRunning diff_rhyth_AD_severity()")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}
  ##### read in ROSMAP clin ####
  rosmap_clin = read_csv(rosmap_clin_path, show_col_types = FALSE)
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

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist

  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    cyc_pred_merged = merge(cyc_pred, rosmap_clin, by.x = "ID", by.y = "projid", y.keep = F)
    preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, batch) %>% arrange(Phase)

  }else{
    cyc_pred_merged = merge(cyc_pred, rosmap_clin, by.x = "ID", by.y = "projid", y.keep = F)
    preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin) %>% arrange(Phase)
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% seedlist), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% seedlist), 1]))  #add the gene names to the columns of gene1

  cog = as.factor(preds$cogdx[match(rownames(gene1), preds$ID)])      # cogdx score 4 or 5
  cerad = as.factor(preds$ceradsc_bin[match(rownames(gene1), preds$ID)])
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  if (useBatch){b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) }

  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I_local_cog = cog
    I_local_cerad = cerad
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I_local_cog = cog[-rm_NA]
        I_local_cerad = cerad[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      if(useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + b1)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1)
      }else{
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + 0)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + 0)
      }

      anova_results = anova(partial_model, full_model)
      p_cog = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]

      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I_local_cog5:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I_local_cog5:cos(times1)"]] + cos_coeff
      acrophase_cog4 = atan2(sin_coeff, cos_coeff) %% (2*pi)
      amplitude_cog4= sqrt((sin_coeff^2) + (cos_coeff^2))
      acrophase_cog5 = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_cog5= sqrt((sin_coeff2^2) + (cos_coeff2^2))

      ####### ceradsc_binned #########
      if(useBatch){
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + b1)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + b1)
      }else{
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + 0)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + 0)
      }
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

mesor_differences = function(cyc_pred, tmm, DR_genes, pb = NULL, useBatch = F, percentile = 0.025){ ##
  cat("\nRunning Mesor_differences()")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist

  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase) %>% arrange(Phase)
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% DR_genes), -1] # "gene" is tmm with only seedlist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% DR_genes), 1]))  #add the gene names to the columns of gene1

  if(useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) # batch
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])  # condition factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])


  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
    if(useBatch){
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + b)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + I + b)
    }else{
      partial_model = lm(gexp1 ~ sin(times) + cos(times))
      full_model = lm(gexp1 ~ sin(times) + cos(times) + I)

    }

    anova_results = anova(partial_model, full_model)
    wilcox_test = wilcox.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    p_wilcox = wilcox_test$p.value

    t_test = t.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    p_ttest = t_test$p.value

    p_mesor = anova_results$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    if(useBatch){
      #When you have multiple batches, which batch do you use as the mesor? We take the weighted average
      # here avg_cond_0_mesor is ((#_cond0_b0 * intercept) + num_cond0_b1 * (intercept + b1_offset)  ) / num_Cond_0
      avg_cond_0_mesor = ( (full_model[["coefficients"]][["(Intercept)"]] * length(which(b == levels(b)[1] & I == levels(I)[1])) ) +
                            ((full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["bcond_1"]]) * length(which(b == levels(b)[2] & I == levels(I)[1])) ) ) /
        length(which(I == levels(I)[1]))

      # here avg_cond_1_mesor is (#_cond1_b0 * (intercept + cond1_offset) + #_cond1_b1 * (intercept + cond1_offset + b1_offset)  ) / num_Cond_1
      avg_cond_1_mesor = ((full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["Icond_1"]] ) * length(which(b == levels(b)[1] & I == levels(I)[2])) +
                            (full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["bcond_1"]] + full_model[["coefficients"]][["Icond_1"]] ) * length(which(b == levels(b)[2] & I == levels(I)[2])) ) /
        length(which(I == levels(I)[2]))
      mesor_AD = avg_cond_1_mesor
      mesor_CTL = avg_cond_0_mesor

      # mesor_prcnt_change = avg_cond_1_mesor / avg_cond_0_mesor
    }else{
      mesor_CTL = full_model[["coefficients"]][["(Intercept)"]]
      mesor_AD = full_model[["coefficients"]][["(Intercept)"]] + full_model[["coefficients"]][["Icond_1"]]

    }
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

    info = cbind( Gene_Symbols, p_mesor, p_wilcox, p_ttest, mesor_CTL, mesor_AD)
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

run_cycling_and_dr_analysis = function(order_path, tmm_path, isCyclingSigCutoff = 0.05, useBatch = F, percentile = 0.025){
  tmm = read_csv(tmm_path, show_col_types = FALSE)      #read expression data, unordered
  colnames(tmm)[1] = "gene_names" #set first column name bc sometimes they are different

  #grab these for comparison to what cyclops finds:
  cyc_pred_file = list.files(path = paste0(order_path, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(order_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)

  #perform is_cycling (method 1) nested regression on CTL data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_CTL = is_cycling(cyc_pred, tmm, cond_subset = "cond_0", useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_CTL$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  cycling_in_CTL = cbind(Ensembl, cycling_in_CTL)
  #record strong cyclers in CTL, for several different amplitude cutoffs
  strong_cyclers_CTL_AR25 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.25 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_CTL_AR33 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.33 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_CTL_AR1 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.1 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  #perform is_cycling (method 1) nested regression on AD data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_AD = is_cycling(cyc_pred, tmm, cond_subset = "cond_1", useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_AD$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  cycling_in_AD = cbind(Ensembl, cycling_in_AD)
  #record strong cyclers in AD, for several different amplitude cutoffs
  strong_cyclers_AD_AR25 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.25 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_AD_AR33 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.33 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_AD_AR1 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.1 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  #perform is_cycling (method 2) nested regression on all data ( all data tested together, same to compareRhythms)
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_either_cond = is_cycling_method2(cyc_pred, tmm, useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_either_cond$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  cycling_in_either_cond = cbind(Ensembl, cycling_in_either_cond)
  #record strong cyclers in either condition, for several different amplitude cutoffs
  strong_cyclers_method2_AR25 = dplyr::filter(cycling_in_either_cond, ( (as.numeric(amp_ratio_CTL) >=0.25) | (as.numeric(amp_ratio_AD) >=0.25)) & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_method2_AR1 = dplyr::filter(cycling_in_either_cond, ( (as.numeric(amp_ratio_CTL) >=0.1) | (as.numeric(amp_ratio_AD) >=0.1)) & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  # We only test for diff rhythmicity if a gene cycles in AD OR CTL. Here I create those unions
  seedlist_AR25 = union(strong_cyclers_AD_AR25$Gene_Symbols, strong_cyclers_CTL_AR25$Gene_Symbols)
  seedlist_AR33 = union(strong_cyclers_AD_AR33$Gene_Symbols, strong_cyclers_CTL_AR33$Gene_Symbols)
  seedlist_AR1 = union(strong_cyclers_AD_AR1$Gene_Symbols, strong_cyclers_CTL_AR1$Gene_Symbols)
  seedlist_method2_AR25 = strong_cyclers_method2_AR25$Gene_Symbols
  seedlist_method2_AR1 = strong_cyclers_method2_AR1$Gene_Symbols

  DR_genelist_list = list(seedlist_AR25, seedlist_AR33, seedlist_AR1, seedlist_method2_AR25, seedlist_method2_AR1)

  ##### mesor differences ######
  gene_list_mesor =  unlist(unname(tmm[!grepl("_D", unlist(tmm[,1])), 1])) # TEST ALL genes for Mesor diff (not just cyclers)
  pb <- progress_bar$new(total = length(gene_list_mesor))
  differential_mesor = mesor_differences(cyc_pred, tmm, gene_list_mesor,useBatch = useBatch, pb = pb, percentile = percentile)
  Ensembl = Ensembl_dict$ENSEMBL[match(differential_mesor$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  differential_mesor = cbind(Ensembl, differential_mesor)

  ####### differential rhtyhms with continuous cerad covs####
  diff_rhythms_AD_severity = diff_rhyth_AD_severity(cyc_pred, tmm,
                                                    unname(unlist(strong_cyclers_AD_AR25$Gene_Symbols)),
                                                    rosmap_clin_path = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/Meta_data/cleaned_rosmap_meta_cogdxConds.csv",
                                                    percentile = percentile, useBatch = useBatch)

  ####### differential rhythms #####
  DR_results = list()
  for(genelist in DR_genelist_list){
    pb <- progress_bar$new(total = length(genelist))
    diff_rhythms_results = diff_rhyth(cyc_pred, tmm, genelist, useBatch = useBatch, pb = pb, percentile = percentile)
    Ensembl = Ensembl_dict$ENSEMBL[match(diff_rhythms_results$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
    diff_rhythms_results = cbind(Ensembl, diff_rhythms_results)
    DR_results = c(DR_results, list(diff_rhythms_results))
  }

  diff_rhythms25 = DR_results[[1]]
  diff_rhythms33 = DR_results[[2]]
  diff_rhythms1 = DR_results[[3]]
  diff_rhythms_mthd2_AR25 = DR_results[[4]]
  diff_rhythms_mthd2_AR1 = DR_results[[5]]

  #Create list of strong cyclers (AR 0.25 or 0.33 and BHQ < BHQcutoff) in CTL subjects
  CTL_cyclers_AR25BHQCutoff = dplyr::select(strong_cyclers_CTL_AR25, Ensembl , Gene_Symbols )
  CTL_cyclers_AR33BHQCutoff = dplyr::select(strong_cyclers_CTL_AR33, Ensembl , Gene_Symbols )

  #Create list of strong cyclers (AR 0.25 or 0.33 and BHQ < BHQcutoff) in AD subjects
  AD_cyclers_AR25BHQCutoff = dplyr::select(strong_cyclers_AD_AR25, Ensembl , Gene_Symbols )
  AD_cyclers_AR33BHQCutoff = dplyr::select(strong_cyclers_AD_AR33, Ensembl , Gene_Symbols )

  #Create list of strong cyclers (AR 0.25 or 0.33 and BHQ < BHQcutoff) in Either (method 2)
  mthd2_cyclers_AR25BHQCutoff = dplyr::select(strong_cyclers_method2_AR25, Ensembl, Gene_Symbols)
  mthd2_cyclers_AR1BHQCutoff = dplyr::select(strong_cyclers_method2_AR1, Ensembl, Gene_Symbols)

  # All genes expressed in CTL and AD:
  EnrichR_background = dplyr::select(cycling_in_CTL, Ensembl, Gene_Symbols)

  #Create list of strong DR genes (AR 0.33 or 0.25 or 0.1, and BHQ< 0.2 or BHQ< 0.1, respectively)
  DR_cyclers_AR33_DRBHQ2 = dplyr::filter(diff_rhythms33,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_AR25_DRBHQ2 = dplyr::filter(diff_rhythms25,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_AR1_DRBHQ2 = dplyr::filter(diff_rhythms1,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_mthd2_AR25_DRBHQ2 = dplyr::filter(diff_rhythms_mthd2_AR25,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))
  DR_cyclers_mthd2_AR1_DRBHQ2 = dplyr::filter(diff_rhythms_mthd2_AR1,  as.numeric(BHQ) < 0.2) %>% arrange(as.numeric(BHQ))

  #Create list of lost cycling DR genes:
  DR_lostAmpAD_AR33BHQ2 = filter(diff_rhythms33, BHQ < 0.2, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_lostAmpAD_AR25BHQ2 = filter(diff_rhythms25, BHQ < 0.2, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_lostAmpAD_AR1BHQ2 = filter(diff_rhythms1, BHQ < 0.2, Log_AD_CTL_ampRatio < 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_mthd2_lostAmpAD_AR25BHQ2 = filter(diff_rhythms_mthd2_AR25, BHQ < 0.2, Log_AD_CTL_ampRatio < 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_mthd2_lostAmpAD_AR1BHQ2 = filter(diff_rhythms_mthd2_AR1, BHQ < 0.2, Log_AD_CTL_ampRatio < 0) %>% dplyr::select(Ensembl, Gene_Symbols)

  #Create lists of gained cycling DR genes
  DR_gainAmpAD_AR33BHQ2= filter(diff_rhythms33, BHQ < 0.2, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_gainAmpAD_AR25BHQ2= filter(diff_rhythms25, BHQ < 0.2, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_gainAmpAD_AR1BHQ2 = filter(diff_rhythms1, BHQ < 0.2, Log_AD_CTL_ampRatio > 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_mthd2_gainAmpAD_AR25BHQ2 = filter(diff_rhythms_mthd2_AR25, BHQ < 0.2, Log_AD_CTL_ampRatio > 0) %>% dplyr::select(Ensembl, Gene_Symbols)
  DR_mthd2_gainAmpAD_AR1BHQ2 = filter(diff_rhythms_mthd2_AR1, BHQ < 0.2, Log_AD_CTL_ampRatio > 0) %>% dplyr::select(Ensembl, Gene_Symbols)


  if (!(dir.exists(paste(order_path, "downstream_output", sep = '/')))){
    dir.create(paste(order_path, "downstream_output", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_results", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_files", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "PSEA_files", sep = '/'))

  }
  #Create string for isCyclingSigCutoff.  E.g. 0.05 -> "05"
  isCyclingSigCutoff_str = str_extract(as.character(isCyclingSigCutoff), "(?<=\\.)\\d+")
  blunting_percentile_str = str_extract(as.character(percentile), "(?<=\\.)\\d+")

  #write out all results of cycling and DR analysis
  write.table(diff_rhythms25, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio25.csv"), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms33, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio33.csv"), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms1, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio1.csv"), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms_mthd2_AR25, paste0(order_path, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio25.csv"), sep = ',', row.names = F, col.names = T)
  write.table(diff_rhythms_mthd2_AR1, paste0(order_path, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio1.csv"), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_CTL, paste(order_path, "downstream_output","cosinor_results_CTL.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_AD, paste(order_path, "downstream_output","cosinor_results_AD.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_either_cond, paste(order_path, "downstream_output","cosinor_results_method2_cyclingInEither.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #gene lists for enrichR
  write.table(CTL_cyclers_AR25BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR25BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(CTL_cyclers_AR33BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR33BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR25BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR25BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR33BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR33BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(mthd2_cyclers_AR25BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_CTL_mthd2_cyclers_AR25BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)

  #background genes for enrichR
  write.table(EnrichR_background, paste(order_path, "downstream_output", "enrichR_files","EnrichR_background.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #DR genes for enrichR
  write.table(DR_cyclers_AR33_DRBHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_AR25_DRBHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_AR1_DRBHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_mthd2_AR25_DRBHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_mthd2_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_cyclers_mthd2_AR1_DRBHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_mthd2_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)

  write.table(DR_lostAmpAD_AR33BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_lostAmpAD_AR25BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_lostAmpAD_AR1BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_mthd2_lostAmpAD_AR25BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_mthd2_lostAmpAD_AR1BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)

  write.table(DR_gainAmpAD_AR33BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_gainAmpAD_AR25BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_gainAmpAD_AR1BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_mthd2_gainAmpAD_AR25BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR25_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)
  write.table(DR_mthd2_gainAmpAD_AR1BHQ2, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ2.csv"), sep = ',', row.names = F, col.names = T)

  #Mesor differences
  write.table(differential_mesor, paste(order_path, "downstream_output", "differential_mesor_all_genes.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  sig_diff_mesor = filter(differential_mesor, as.numeric(BHQ) < 0.05 ) %>% dplyr::select( Ensembl, Gene_Symbols)
  write.table(sig_diff_mesor, paste(order_path, "downstream_output","enrichR_files","diff_mesor_all_genes_BHQ05.csv", sep = "/"), sep = ',', row.names = F, col.names = T)

  #Continuous AD differences
  write.table(diff_rhythms_AD_severity, paste(order_path, "downstream_output", "diff_rhythms_AD_severity_AR25.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  strong_cogdx_diffs = filter(diff_rhythms_AD_severity, BHQ_cogdx< 0.1) %>% dplyr::select(Gene_Symbols)
  Ensembl = Ensembl_dict$ENSEMBL[match(strong_cogdx_diffs$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  strong_cogdx_diffs = cbind(Ensembl, strong_cogdx_diffs) %>% dplyr::select(Ensembl , Gene_Symbols )
  write.table(strong_cogdx_diffs, paste(order_path, "downstream_output", "enrichR_files", "strong_cogdx_diffs_AR25.csv", sep = "/"), sep = ',', row.names = F, col.names = T)

  #create lists of genes for PSEA
  PSEA_CTL_cyclers_AR25BHQCutoff = strong_cyclers_CTL_AR25 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_CTL_cyclers_AR25BHQCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_CTL_cyclers_AR25BHQ", isCyclingSigCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  PSEA_AD_cyclers_AR25BHQCutoff = strong_cyclers_AD_AR25 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_AD_cyclers_AR25BHQCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_AD_cyclers_AR25BHQ", isCyclingSigCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  PSEA_DR_AR25BHQ2_acrodiffs = DR_cyclers_AR25_DRBHQ2 %>% mutate(acro_diff = (as.numeric(acrophase_AD) - as.numeric(acrophase_CTL))*12/pi ) %>%
    dplyr::select(Gene_Symbols, acro_diff)
  write.table(PSEA_DR_AR25BHQ2_acrodiffs, paste0(order_path, "/downstream_output/PSEA_files/PSEA_DR_AR25BHQ2_acrodiffs.txt"), sep = '\t', row.names = F, col.names = F, quote = F)

  #write out nice summary of cycling and DR genes
  summary = data.frame(List = c("TMM name","Using batch in Regression", "isCyclingBHQCutoff", "Blunting_Percentile" ,paste0("CTL_cyclers_AR25BHQ", isCyclingSigCutoff_str), paste0("CTL_cyclers_AR33BHQ", isCyclingSigCutoff_str),
      paste0("AD_cyclers_AR25BHQ", isCyclingSigCutoff_str), paste0("AD_cyclers_AR33BHQ", isCyclingSigCutoff_str), paste0("AD_CTL_mthd2_cyclers_AR25BHQ", isCyclingSigCutoff_str),
      "DR_cyclers_AR1BHQ2", "DR_cyclers_AR25BHQ2", "DR_cyclers_AR33BHQ2", "DR_cyclers_mthd2_AR1_DRBHQ2", "DR_cyclers_mthd2_AR25_DRBHQ2",
      "DR_lostAmpAD_AR1BHQ2", "DR_lostAmpAD_AR25BHQ2", "DR_lostAmpAD_AR33BHQ2", "DR_mthd2_lostAmpAD_AR1BHQ2", "DR_mthd2_lostAmpAD_AR25BHQ2",
      "DR_gainAmpAD_AR1BHQ2", "DR_gainAmpAD_AR25BHQ2", "DR_gainAmpAD_AR33BHQ2", "DR_mthd2_gainAmpAD_AR1BHQ2", "DR_mthd2_gainAmpAD_AR25BHQ2"),
      Num_genes = c(tmm_path, useBatch, isCyclingSigCutoff, percentile, dim(CTL_cyclers_AR25BHQCutoff)[1], dim(CTL_cyclers_AR33BHQCutoff)[1],
      dim(AD_cyclers_AR25BHQCutoff)[1],dim(AD_cyclers_AR33BHQCutoff)[1],
      dim(mthd2_cyclers_AR25BHQCutoff)[1],
      dim(DR_cyclers_AR1_DRBHQ2)[1],dim(DR_cyclers_AR25_DRBHQ2)[1], dim(DR_cyclers_AR33_DRBHQ2)[1], dim(DR_cyclers_mthd2_AR1_DRBHQ2)[1], dim(DR_cyclers_mthd2_AR25_DRBHQ2)[1],
      dim(DR_lostAmpAD_AR1BHQ2)[1], dim(DR_lostAmpAD_AR25BHQ2)[1],dim(DR_lostAmpAD_AR33BHQ2)[1], dim(DR_mthd2_lostAmpAD_AR1BHQ2)[1], dim(DR_mthd2_lostAmpAD_AR25BHQ2)[1],
      dim(DR_gainAmpAD_AR1BHQ2)[1], dim(DR_gainAmpAD_AR25BHQ2)[1], dim(DR_gainAmpAD_AR33BHQ2)[1], dim(DR_mthd2_gainAmpAD_AR1BHQ2)[1], dim(DR_mthd2_gainAmpAD_AR25BHQ2)[1])
  )
  write.table(summary, paste(order_path, "downstream_output","cosinor_DR_summary.csv", sep = '/'), sep = ',', row.names = F, col.names = T)


}
