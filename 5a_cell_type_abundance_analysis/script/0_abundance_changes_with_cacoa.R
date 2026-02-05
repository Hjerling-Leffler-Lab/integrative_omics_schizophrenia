#author: Lisa Bast
#date: 13.09.2022
#about: cacoa for snRNAseq SCZ CTRL data (https://github.com/kharchenkolab/cacoa) to perform comparitive analysis of case vs. control
#       i) analysis of compositional shifts (cluster-based/ using CT annotation)
#       ii) analysis of transcriptional shifts (cluster-free/ not using CT annotation)
#       Fig. 1F, S5A
#version: 0.0.3

library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(cacoa)
library(trqwe)
library(coda.base)
library(reshape2)
library(ggpubr)
library(ggbreak)
library(cowplot)
library(SeuratObject)
library(Seurat)
library(openxlsx)
library(loomR)
library(varhandle) #unfactor()
library(ape)


#settings:
opt_run_per = "celltype" # "all" #"class" #
#opt_inhibitory_only <- FALSE#TRUE #if TRUE runs for all inhibitory subtypes on the 37 cluster level
data_subset_modes = c("inhibitory_only","excitatory_only","non_neuronal_only","neurons_only")#,"all_cells")

opt_species='human'
comp_name = 'longleaf'#"PowerWS" #

if (comp_name == 'longleaf'){
  opt_run_data_subset_modes_in_parallel = TRUE
} else{
  opt_run_data_subset_modes_in_parallel = FALSE
}

n_cores <- 35

code_path = getwd()
source("utils.R")

setwd("../")
project_path = getwd()
setwd("../")
main_path = getwd()

n_cluster=15 #37

if (opt_run_data_subset_modes_in_parallel){
  args <- as.integer(commandArgs(trailingOnly=TRUE)) # ids referring to data subset mode for running them in parallel
} else{
  args <- 1 # for first option
}

data_subset_mode <- data_subset_modes[args]

## define results path
results_path <- paste0(project_path,"/results/Cell_type_abundance/",data_subset_mode,"/")

#make sure results path exists
if (file.exists(results_path)==FALSE) {
    dir.create(results_path, recursive=TRUE)
}


#Input Loom file(): Ct annotated and clustered cell type classes, downsampled to 20% of cells for testing, otherwise all cells
data_path<-paste0(main_path,'/4_data_integration_and_cell_type_annotation/output/')
if (comp_name == 'longleaf'){
  file_name_add_str <- ""
} else{
  file_name_add_str <- "_subsampled_to_10_percent_cells"
}

#get color_palettes:
palettes <- get_color_palettes(data_path)

#specify loom file 
if (data_subset_mode == "inhibitory_only"){
  filename <- paste0("Samples_Inhibitory_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered",file_name_add_str,".loom")
}else if (data_subset_mode == "excitatory_only"){
  filename <- paste0("Samples_Excitatory_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered",file_name_add_str,".loom")
}else if(data_subset_mode == "non_neuronal_only"){
  filename <- paste0("Samples_NonNeuronal_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered",file_name_add_str,".loom")
}else {
  filename <- paste0("Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered",file_name_add_str,".loom")
}

#load raw expression count matrix
data <- load_data_for_cacoa(data_path,filename)
print("loading data works")

cao <- build_cacoa_object(data,n_cluster,n_cores,data_subset_mode,palettes,opt_run_per)

cluster_based_analysis(cao, results_path)


