#####
# Tues Sept 5 2023
# Script for taking a df containing TF_names (from pscan but also enrichR results),
# and checking if they are found in DE, DE, DM, cycling files. 

TF_filename = paste0(path_to_cyclops_ordering, "diff_rhythms/pscan/pscan_fDR_cyclers_AR1BHQ1_pscan.csv") #TODO remove this line, add as an arg
#open DR genes
TF_names = read.csv(TF_filename)

DR_AR1 = read.csv(paste0(path_to_cyclops_ordering, "diff_rhythms/diff_rhythms_AmpRatio1.csv"))
