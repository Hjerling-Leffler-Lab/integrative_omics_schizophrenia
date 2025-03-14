#author: Lisa Bast
#date: 07.12.2020, last update 05.07.2024
#version: 0.1.2
#about: automatic cell type annotation with scmap
#       can optionally perform a test on a subset of the reference data set (train and test split up)


library(SingleCellExperiment)
library(scmap)
#library(Seurat) # required for as.sparse method
library(iotools) # requires for read.csv.raw efficiently read csv files
library(pryr) # to analyse memory
library(pheatmap)
library(irr)#kappa calculation
library(doParallel)
library(foreach)
library(loomR)
library(stringr)

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path_project = getwd()

## settings:
opt_server =  TRUE #FALSE# 

alignment_method = 'cellranger'
opt_load_from_loom = TRUE #FALSE
opt_save_as_loom = TRUE# FALSE
opt_doublet_version = '_'
#should high drop out rate genes be specified?
opt_determine_HDR_genes=FALSE

ref_data = "Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq"
ref_data_folder = "ABM_MCA/"
ref_data_subfolder_str="red_to_"

#if (opt_server==TRUE){
#  args <- commandArgs(trailingOnly=TRUE)
#}else{
#  args <- "3"
#}
#n <- as.integer(args)

data_sets <- c("Bast_et_al","Velmeshev_et_al","ROSMAP_Inhibitory_Neurons_1","ROSMAP_Excitatory_Neurons_1","ROSMAP_Inhibitory_Neurons_2","ROSMAP_Excitatory_Neurons_2")
for (n in seq(2,6)){
  
	opt_annotate <- data_sets[n]
	print(opt_annotate)
	if (opt_annotate=="Bast_et_al"){
	  n_CT = 3#4 #4 or 51 or 76 or 120
	} else if (opt_annotate=="Velmeshev_et_al"){
	  n_CT = 51# or 76 or 120
	} else{ # ROSMAP
	  n_CT = 51
	}

	# setup parallel backend to use many processors
	#cl <- makeCluster(max(1,detectCores()-1)) #not to overload your computer
	#registerDoParallel(cl)

	opt_scmap <- 'cell2cluster'#'cluster' #
	n_cells_query <- c(2500) # number of cells in query files
	n_genes_ref <- c(5000) #c(500) #number of genes selected for reference data set
	red_method <- c('HVG_with_cell_ranger') #method used for reducing number of genes in reference data set
	ref_data_subfolder_str <- paste0(ref_data_subfolder_str,as.character(n_genes_ref),"_",red_method,"/")

	## 1) specify paths and load function to read data
	path_code <- paste0(main_path_project,'/script/')

	path_results <- paste0(main_path_project,'/output/CT_anno_scmap/Filtered_',ref_data_folder,alignment_method,'/')

	path_results_figures <- paste0(path_results,n_CT,"_CTs/",'figures',opt_doublet_version,'/')
	path_results_files <- paste0(path_results,n_CT,"_CTs/",'output_files',opt_doublet_version,'/')

	## 2) specify file names
	if (opt_annotate=="Bast_et_al"){
	  path_query_data <- paste0(main_path_project,'/data/')
	  query_file = paste0(path_query_data,"Samples_cellranger_pagoda_TH_and_D_adj_filtered_conos_cluster_based_Endothelial_and_mural_cells.loom")
	} else if (opt_annotate=="Velmeshev_et_al"){
	  path_query_data <- paste0(main_path_project,'/data/') #/cfs/klemming/projects/supr/snic2020-6-62/SCZ_human_sc_PFC/results/MultiNeuronChat/Velmeshev_2019/data
	  query_file = paste0(path_query_data,"Velmeshev_2019_matrix.loom")
	} else if (opt_annotate=="ROSMAP_Excitatory_Neurons_1"){
	  path_query_data <- paste0(main_path_project,'/data/') 
	  query_file = paste0(path_query_data,"joined_matrix_subsampled_Excitatory_Neurons_1.loom")
	} else if (opt_annotate=="ROSMAP_Excitatory_Neurons_2"){
	  path_query_data <- paste0(main_path_project,'/data/') 
	  query_file = paste0(path_query_data,"joined_matrix_subsampled_Excitatory_Neurons_2.loom")
	} else if (opt_annotate=="ROSMAP_Inhibitory_Neurons_1"){
	  path_query_data <- paste0(main_path_project,'/data/') 
	  query_file = paste0(path_query_data,"joined_matrix_subsampled_Inhibitory_Neurons_1.loom")
	} else if (opt_annotate=="ROSMAP_Inhibitory_Neurons_2"){
	  path_query_data <- paste0(main_path_project,'/data/') 
	  query_file = paste0(path_query_data,"joined_matrix_subsampled_Inhibitory_Neurons_2.loom")
	}

	ref_file_symbol_name <- paste0('9606_symbol_',n_CT,'.csv')
	ref_file_map_name <- paste0('9606_map_',n_CT,'.csv')
	subpath_ref_data <- paste0(main_path_project,'/data/reference_data_sets/',ref_data_folder,sub("_2","",sub("_1","",opt_annotate)),"_",ref_data_subfolder_str)
	paste("_",paste(substr(query_file,1,3),substr(query_file,nchar(query_file)-5,nchar(query_file)-4),sep=""),sep="")

	print(query_file)
	print(subpath_ref_data)

	## 3) annotate cell types
	annotate_cell_types_with_scmap(red_method, 
								   n_genes_ref, 
								   subpath_ref_data, 
								   ref_file_symbol_name,
								   ref_data_folder,
								   opt_load_from_loom,
								   opt_save_as_loom,
								   opt_doublet_version,
								   path_results_figures,
								   n_CT,
								   query_file,
								   path_results_files)
	  
	#stopCluster(cl)
}

