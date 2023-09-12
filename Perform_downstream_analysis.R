library(reticulate) #use to find conda env for python
library(tidyverse)

#get path to where this file resides
path = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path) #setwd to path of this file
source("./differential_rhyth_single_cov_for_AnalysisPipeline.R")
source("./create_rnk_files.R")
source("./fgsea.R")
source("./compareRhythms.R")
source("./Mesor_Venn_Diagram.R")
source("./Files_for_Pscan.R")
source("./Make_KEGG_map_files.R")
source("./Check_TF_for_cycling_or_DR.R")
conda_list()
conda_list()[[2]][2] %>% 
  use_condaenv(required = TRUE)

# absolute path to translation dict between gene symbols and ENSEMBL
Ensembl_dict = readxl::read_xlsx("scROSMAP_ENSEMBL_dict.xlsx")
entrez_dict = readxl::read_excel("scROSMAP_ENSEMBL_ENTREZ_dict.xlsx")

#optional, add processors with doParallel pkg
# cl <- makePSOCKcluster(8)
# registerDoParallel(cl)

run_downstream_analysis = function(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.05, deseq_de_filename = NULL){
  # this function handles, cosinor regression (custom and CompareRhtyms),
  # for cycling and differential cycling, EnrichR and fGSEA pathway analysis.
  # input: path_to_cyclops_ordering: string representing absolute path to cyclops ordering
  #.       path_to_tmm_file: string: absolute path to tmm gene expression matrix (csv)
  
  #run custom cycling and differential cycling analysis
  print("***Running cosinor and differential rhythms analysis***")
  run_cycling_and_dr_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = isCyclingBonfCutoff)
  #Run the cycling and DC analysis with CompareRhythms package
  print("***Running CompareRhythms Package***")
  run_compare_rhythms(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = isCyclingBonfCutoff)
  
  #creates a Venn diagram of diff_mesor genes vs diff mean genes from ttest
  print("***Outputing Venn Diagram***")
  create_mesor_venn_diag(path_to_cyclops_ordering, BHQ = 0.05)
  
  #############################################
  # EnrichR on cycling and DR cycling results #
  #############################################
  setwd(paste0(path_to_cyclops_ordering, "diff_rhythms"))
  path_no_space = gsub(" ", "\\ ", path, fixed = TRUE)
  
  gene_lists = list.files("./enrichR_files/", pattern = "^(CTL_|AD_|DR_)")
  sapply(gene_lists, function(x){ 
    system(paste0("python3 ", path_no_space, "/Python_EnrichR_for_AnalysisPipeline.py -g \"enrichR_files/" ,x, "\" -b \"enrichR_files/EnrichR_background.csv\""))})
  
  setwd("enrichR_files/")
  DR_gene_lists = list.files( pattern = "^(DR_)")
  CR_files = grepl("COMPARERHYTHMS", DR_gene_lists)
  DR_backgrounds = paste0("diff_rhythms_AmpRatio", str_replace(str_extract(DR_gene_lists, pattern = "AR\\d+"), "AR", ""), ".csv")
  DR_backgrounds[CR_files] = "diff_rhythms_AmpRatio25_COMPARERHYTHMS.csv"
  
  ##############
  #Diff_rhythms# 
  ##############
  setwd(paste0(path_to_cyclops_ordering, "diff_rhythms"))
  mapply(function(x, y){ 
    system(paste0("python3 ", path_no_space, "/Python_EnrichR_for_AnalysisPipeline.py -g \"enrichR_files/" ,x, "\" -b \"",y , "\""))},
    DR_gene_lists, 
    DR_backgrounds)
  #AD severity
  system(paste0("python3 ", path_no_space, "/Python_EnrichR_for_AnalysisPipeline.py -g \"enrichR_files/strong_cogdx_diffs_AR25.csv\" -b \"enrichR_files/AD_cyclers_AR25Bonf05.csv\""))
  
  ##########
  # fGSEA  #
  ##########
  setwd(path_to_cyclops_ordering)
  if (!(dir.exists("diff_rhythms/fGSEA"))){
    dir.create("diff_rhythms/fGSEA")
    dir.create(paste("diff_rhythms", "fGSEA", "rnk_files", sep = "/"))
    dir.create(paste("diff_rhythms" , "fGSEA", "fGSEA_results", sep = "/"))
  }
  
  
  #create rnk files
  write_rnks(path_to_cyclops_ordering)
  
  #path of parent folder
  setwd(path)
  
  #Pathways, downloaded from MsigDB
  pathways <- c(gmtPathways("./MsigDB_gmts_for_GSEA/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
                , gmtPathways("./MsigDB_gmts_for_GSEA/h.all.v2023.1.Hs.symbols.gmt"))
  gene_remapping_dict = read.delim("./MsigDB_gmts_for_GSEA/Human_Gene_Symbol_with_Remapping_MSigDB.v2023.1.Hs.chip", sep = '\t')
  setwd(paste(path_to_cyclops_ordering,"diff_rhythms", "fGSEA", "rnk_files", sep = "/"))
  rnk_files = list.files(pattern = ".rnk")
  
  
  if(!(dir.exists(paste(path_to_cyclops_ordering, "diff_rhythms","fGSEA","fGSEA_results", sep = "/")))){
    dir.create( paste(path_to_cyclops_ordering,"diff_rhythms", "fGSEA","fGSEA_results", sep = "/") )
    dir.create( paste(path_to_cyclops_ordering,"diff_rhythms", "fGSEA","fGSEA_results", "plots",sep = "/"))
  }
  
  run_fgsea(rnk_files, gene_remapping_dict, pathways)
  
  
  ########
  # PSEA #
  ########
  setwd(path)
  #Pathways, downloaded from MsigDB, absolute paths
  psea_pathways <- c(paste0(path, "/MsigDB_gmts_for_GSEA/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
                , paste0(path, "/MsigDB_gmts_for_GSEA/h.all.v2023.1.Hs.symbols.gmt"))
  psea_pathways = gsub(" ", "\\ ", psea_pathways, fixed = TRUE) #make sure no spaces when calling cmd
  
  # tab delimed txt files 
  setwd(paste(path_to_cyclops_ordering,"diff_rhythms", "PSEA_files",sep = "/"))
  psea_files = list.files(pattern = "*.txt")
  
  path_to_cyclops_ordering_no_space = gsub(" ", "\\ ", path_to_cyclops_ordering, fixed = TRUE)

  arg_df = expand.grid(psea_files, psea_pathways)%>% 
                mutate(out_dirs = paste(str_extract(string = as.character(Var1), pattern= ".+(?=\\.txt)"), 
                                     basename(as.character(Var2)) %>% str_extract(pattern= ".+(?=\\.gmt)")%>% 
                gsub(pattern = "\\.", replacement = ""), sep = "_") )
  apply(arg_df, 1,
        function(a){
          system(paste0("java -jar ", path_no_space, "/PSEA_cmd.jar ", paste0(path_to_cyclops_ordering_no_space, "/diff_rhythms/PSEA_files/", a[[1]]), " ", a[[2]]," ", paste0(path_to_cyclops_ordering_no_space, "/diff_rhythms/PSEA_files/", a[[3]]), " 5 10000 pdf"))
          })
  #########
  # Pscan #
  #########
  print("Starting Pscan Analysis")
  if(!(dir.exists(paste0(path_to_cyclops_ordering, "diff_rhythms/pscan")))){
    dir.create(paste0(path_to_cyclops_ordering, "diff_rhythms/pscan") )
  }
  setwd(paste0(path_to_cyclops_ordering, "diff_rhythms"))
  
  gene_lists = list.files("./enrichR_files/", pattern = "^(CTL_|AD_|DR_)") 
  lapply(gene_lists, write_pscan_input)
  path_no_space = gsub(" ", "\\ ", path, fixed = TRUE)
  
  setwd("pscan/")
  pscan_files = list.files("./pscan_files/", pattern = ".csv$") 
  
  sapply(pscan_files, function(x){ 
    system(paste0("python3 ", path_no_space, "/pscan_requests.py --file pscan_files/" ,x))})
  
  if(!is.null(deseq_de_filename)){
    setwd("pscan_results")
    pscan_result_files = list.files(pattern = ".csv$")
    print("Searching cycling and DR results for pscan TF")
    sapply(pscan_result_files, augment_tf_file, deseq_de_filename)
    setwd(paste0(path_to_cyclops_ordering, "diff_rhythms/enrichR_results"))
    directories = list.files()
    enrichR_TF_paths <- lapply(directories, function(dir) list.files(dir, full.names = TRUE, recursive = TRUE, pattern = ".*TRANSFACandJASPARPWMs_BCKGRND.*\\.csv|.*TranscriptionFactorPPIs_BCKGRND.*\\.csv"))
    enrichR_TF_paths = flatten(enrichR_TF_paths) %>% paste0("./", .)  
    sapply(enrichR_TF_paths, augment_tf_file, deseq_de_filename)
    
    }
  
  ##############################
  # Write out KEGG Image files #
  ##############################
  setwd(path_to_cyclops_ordering)
  if (!(dir.exists("diff_rhythms/KEGG_map_diagrams"))){
    dir.create("diff_rhythms/KEGG_map_diagrams")
  }
  write_kegg_map_files("diff_rhythms/diff_rhythms_AmpRatio1.csv","Log_AD_CTL_ampRatio", entrez_dict)
}

#Microglia
path_to_cyclops_ordering = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/Microglia/Mglia_CellsFiltered10Percent_OrderedbyExcNeuronsCellsFilteredErikChenZhangMinCV14CondCovs3EG_Jun12Redo/"
path_to_tmm_file = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/Mglia_cellsFiltered1count10percentCells.csv"
deseq_de_filename = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Microglia_CellsFiltered10Percent_cogdx_DE_DEseq2.csv"
run_downstream_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.1)

#Astrocytes
path_to_cyclops_ordering = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/Astrocytes/Astro_OrderingFromExcitatoryNeuronsErikChenZhangCondCovs3EG_jun12Redo/"
path_to_tmm_file = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/Astrocyte_cellsFiltered1count10percentCells.csv"
deseq_de_filename = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/Astrocyte_CellsFiltered10Percent_cogdx_DE_DEseq2.csv"
run_downstream_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.1)

#Excitatory Neurons
path_to_cyclops_ordering = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/ExcitatoryNeurons/Exc_Neurons_CellFiltered10Percent_ErikChenZhang_condCovs_3EG_Jun12Redo/"
path_to_tmm_file = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/ExcNeurons_cellsFiltered1count10percentCells.csv"
deseq_de_filename = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/ExcitatoryNeurons_CellsFiltered10Percent_cogdx_DE_DEseq2.csv"
run_downstream_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.05, deseq_de_filename)

#Inhibitory Neurons
path_to_cyclops_ordering = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/InhibitoryNeurons/InhNeurons_CellsFiltered10Percent_OrderedbyExcNeuronsCellsFilteredErikChenZhangMinCV14CondCovs3EG/"
path_to_tmm_file = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/InhNeurons_cellsFiltered1count10percentCells.csv"
deseq_de_filename = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/InhibitoryNeuron_CellsFiltered10Percent_cogdx_DE_DEseq2.csv"
run_downstream_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.05)

#Pseudobulk All Celltypes
path_to_cyclops_ordering = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/training_output/scROSMAP/cogdx_controls/wAD/Pseudobulk/Pseudobulk_10countsIn10Prcnt_OrderedByExcNeuronsCellsFilteredErikChenZhangcondCovs3EG/"
path_to_tmm_file = "~/Box Sync/Henry_stuff/AD_project/human_data/Cyclops_folders/tmms/scROSMAP/cogdx_controls/PseudoBulk_scrosmapALL_cogdxControls_condCovs_filtered10counts10prcntSubs.csv"
#deseq_de_filename = "~/Box Sync/Henry_stuff/AD_project/scROSMAP/simple_differential_expr/InhibitoryNeuron_CellsFiltered10Percent_cogdx_DE_DEseq2.csv"
run_downstream_analysis(path_to_cyclops_ordering, path_to_tmm_file, isCyclingBonfCutoff = 0.05)
