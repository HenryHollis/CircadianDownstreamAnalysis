library(devtools)
library(readxl)
library(tidyverse)
library(stringr)
library(compareRhythms)

## I stopped blunting outliers for compareRhythms because 
## logging data is already doing that

# blunt_outliers = function(vec){
#   num = length(which(!is.na(vec)))
#   ord = sort(vec)
#   upper_val = ord[as.integer(num-num*.025)]
#   lower_val = ord[as.integer(num-num*.975)]
#   
#   vec[which(vec > upper_val)] = upper_val
#   vec[which(vec < lower_val)] = lower_val
#   return(vec)
# }

run_compare_rhythms = function(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.05 ){

TMM = read_csv(path_to_tmm_file)
setwd(paste0(path_to_cyclops_ordering, "/Fits"))
fit = list.files(pattern = "Fit_Output")
cyclops_ordering = read_csv(fit)
times = cyclops_ordering$Phase * 12 / pi

emat = TMM[-c(1:3),]
emat = apply(emat[,-1], 2, as.numeric)
# emat = t(apply(emat, 1, blunt_outliers))

rownames(emat) = TMM[-c(1:3), 1] %>% unlist
cond = as.factor(as.character(TMM[1,-1]))
cond = relevel(cond, "cond_0")


exp_design = data.frame(group = cond, time = times)
logged_emat = log2(emat+1) #add one to keep in domain of Log

print(paste("Running CompareRhythms with isCyclingBonfCutoff =", isCyclingBonfCutoff))
CR_cosinor_method = compareRhythms(logged_emat, exp_design, method = 'cosinor',
                                   just_classify = F, compare_fdr = 0.2,
                                   amp_cutoff = 0.737, rhythm_bonf_cutoff = isCyclingBonfCutoff)
#CR_modsel_method = compareRhythms(logged_emat, exp_design, method = 'mod_sel', just_classify = F)
#table(CR_cosinor_method$category)
#table(CR_modsel_method$category)
colnames(CR_cosinor_method)[1] = "Gene_Symbols"
CR_cosinor_method$Gene_Symbols = gsub("\\.", "\\-", CR_cosinor_method$Gene_Symbols)
cycling_in_CTL = filter(CR_cosinor_method, rhythmic_in_cond_0 == T) %>% dplyr::select(Gene_Symbols)
cycling_in_AD = filter(CR_cosinor_method, rhythmic_in_cond_1 == T) %>% dplyr::select(Gene_Symbols)
DR_cyclers = filter(CR_cosinor_method, diff_rhythmic == T) %>% dplyr::select(Gene_Symbols)
gain_cyclers = filter(CR_cosinor_method, diff_rhythmic == T & category == "gain") %>% dplyr::select(Gene_Symbols)
loss_cyclers = filter(CR_cosinor_method, diff_rhythmic == T & category == "loss") %>% dplyr::select(Gene_Symbols)


if (!(dir.exists(paste(path_to_cyclops_ordering, "downstream_output", sep = '/')))){
  dir.create(paste(path_to_cyclops_ordering, "downstream_output", sep = '/'))

}
#Create string for isCyclingBonfCutoff.  E.g. 0.05 -> "05"
isCyclingBonfCutoff_str = str_extract(as.character(isCyclingBonfCutoff), "(?<=\\.)\\d+")

#write out all results of cycling and DR analysis
write.table(cycling_in_CTL, paste0(path_to_cyclops_ordering, "downstream_output/enrichR_files/CTL_cyclers_COMPARERHYTHMS_AR25Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
write.table(cycling_in_AD, paste0(path_to_cyclops_ordering, "downstream_output/enrichR_files/AD_cyclers_COMPARERHYTHMS_AR25Bonf", isCyclingBonfCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
write.table(DR_cyclers, paste(path_to_cyclops_ordering, "downstream_output", "enrichR_files", "DR_cyclers_COMPARERHYTHMS_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
write.table(gain_cyclers, paste(path_to_cyclops_ordering, "downstream_output", "enrichR_files", "DR_gainAmpAD_COMPARERHYTHMS_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
write.table(loss_cyclers, paste(path_to_cyclops_ordering, "downstream_output", "enrichR_files", "DR_lostAmpAD_COMPARERHYTHMS_AR25BHQ2.csv", sep = '/'), sep = ',', row.names = F, col.names = T)


write.table(CR_cosinor_method, paste(path_to_cyclops_ordering, "downstream_output", "diff_rhythms_AmpRatio25_COMPARERHYTHMS.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
summary = data.frame(List = c(paste0("cycling_in_CTL_AR25Bonf", isCyclingBonfCutoff_str), paste0("cycling_in_AD_AR25Bonf", isCyclingBonfCutoff_str),
                              "COMPARERHYTHMS_DR_cyclers_AR25BHQ2", "COMPARERHYTHMS_gain_cyclers_AR25BHQ2", "COMPARERHYTHMS_lose_cyclers_AR25BHQ2"),
                     Num_genes = c( dim(cycling_in_CTL)[1], dim(cycling_in_AD)[1],
                                    dim(DR_cyclers)[1], dim(gain_cyclers)[1], dim(loss_cyclers)[1]))

write.table(summary, paste(path_to_cyclops_ordering, "downstream_output", "cosinor_DR_summary.csv", sep = '/'), sep = ",", 
            append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE) 


}

