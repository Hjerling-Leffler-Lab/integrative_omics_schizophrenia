# runs hdWGCNA (https://github.com/smorabit/hdWGCNA, https://www.biorxiv.org/content/10.1101/2022.09.22.509094v1.full.pdf),
# some functions are slightly adapted to output additional values or change plotting functions, maths are not touched
# for each cell type build network based on snRNAseq case control data set and test for differential expression, perform pathway analysis for each module
# Fig. S19A-C, 5A
# author: Lisa Bast 
# date: 11th April 2023
# version: 0.0.1

# single-cell analysis package
library(Seurat)
library(SeuratDisk)
library(harmony)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(gridExtra)
library(ggforestplot)
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

#pathway mapping:
library(gprofiler2)

# compute marker gene overlaps
library(GeneOverlap)

#time code --> improve efficiency
library(tictoc)

#load functions (should be in same folder as this script):
code_path = getwd()
source("utils_hdWGCNA.R")
source("utils_settings.R")
source("utils_plotting.R")

##### SETTINGS #####
#to do: specify
num_threads = 24

#opt_save_obj <- FALSE
opt_run_gprofiler <- TRUE
opt_plot_dendro_heatmap <- FALSE#TRUE # this is very slow
opt_load_seurat_obj_with_metacells <- FALSE#TRUE#
opt_save_seurat_obj_with_metacells <- FALSE#TRUE#
opt_filter_genes <-"variable"#, "fraction", "all", or "custom"
opt_save_seurat_obj <- FALSE#TRUE to fine tune network plots and differential expression analysis without rerunning everything, file is tooo big, does not work properly!
test_options <- c("LR") #,"wilcox","bimod","roc","t","negbinom","poisson")#,"MAST" package not installed#,"DESeq2" needs unnormalized counts matrix # as in FindMarkers function in seurat, which is called by DME analysis


if (opt_load_seurat_obj_with_metacells == TRUE){
  library(foreach)
  #cell type id is input parameter
  CT_id <- commandArgs(trailingOnly=TRUE) # takes the input value specified in .sh file that runs this function
} else{
  CT_id <- 2
}
#get data settings
opt_data <- get_opt_data_hdWGCNA()
#get WGCNA settings
settings <- get_settings_hdWGCNA(opt_data)
#get paths
setwd("../")
setwd("../")
opt_data$main_path  = getwd()
opt_data <- get_paths_hdWGCNA(opt_data, settings, opt_filter_genes)

#gene signature score computation
if (settings$gene_sign_score_comp == 'UCell'){
  library(UCell)
}

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
if (num_threads>1){
  enableWGCNAThreads(nThreads = num_threads)
}

#compute matecells
L <- get_metacells(opt_load_seurat_obj_with_metacells, opt_data, settings, opt_filter_genes)
seurat_obj <- L[[1]]
celltype_list <- L[[2]]
n_metacells_per_group <- L[[3]]
network_gene_list <- L[[4]]

#get cell type ID and name
i <- as.integer(CT_id)
ct <- as.character(celltype_list[i])
ct_ <- gsub(" ", "_", ct)
ct_ <- gsub("/", "-", ct_)
#print cell type:
print(ct)

#update dendro settings
settings <- update_dendrogram_settings(settings, ct, opt_filter_genes)

# filtering step to remove cell types that have too few or no metacells 
if (any(row.names(n_metacells_per_group) == ct) & all(!n_metacells_per_group[ct,] == 0)) {
  print(paste0("Sufficient number of metacells for both groups. Continuing with network analysis for ", ct))
} else {
  print(paste0("One group has zero metacells. Skipping network analysis for ", ct))
  next
}

if (opt_filter_genes=="custom"){
  #acvtivate a different set of genes but same metacells as before (created based on network genes for any of the cell types)
  seurat_obj <- SetupForWGCNA(
                  seurat_obj,
                  gene_select = opt_filter_genes,
                  gene_list = network_gene_list[[ct]],
                  wgcna_name = ct,
                  metacell_location = "all_CTs_all_network_genes")
} else{
  #acvtivate a different set of genes but same metacells as before (created based on network genes for any of the cell types)
  seurat_obj <- SetupForWGCNA(
                  seurat_obj,
                  gene_select = opt_filter_genes, #How to select genes? Select "variable", "fraction", "all", or "custom".
                  wgcna_name = ct,
                  metacell_location = "all_CTs_all_network_genes")
}

opt_data$results_path_ct_norm <- create_ct_subfolder(opt_data,ct_)
browser()
#calculate network
seurat_obj <- network_analysis(seurat_obj,ct,ct_,opt_data,settings)

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=opt_data$sample_variable
)
seurat_obj <- module_analyis(seurat_obj,opt_data,ct,settings)

if (opt_save_seurat_obj==TRUE){
  #save seurat object:
  save_seurat_obj(seurat_obj,opt_data,ct_)
}
#visualization:
#1) the network
seurat_obj <- visualize_network(seurat_obj, opt_data$results_path_ct_norm)

#2) the dendrogram
if (opt_plot_dendro_heatmap) {
  plot_dendro_heatmap(seurat_obj, opt_data$results_path_ct_norm, ct)
}

#pathway analysis
if (opt_run_gprofiler){
  run_gprofiler_for_network_modules(seurat_obj,ct_,opt_data)
}

#identify which modules are differentially expressed between the groups
for (opt_test in test_options){
  print(opt_test)
  DME_list <- list()
  tryCatch(
    {
      DME_list <- DME_analysis(seurat_obj, ct, ct_, opt_data, settings, DME_list, opt_test)
  	  DMEs <- do.call(rbind, DME_list)
  	  DMEs$level <- factor(as.character(DMEs$level), levels = celltype_list)
  	  
  	  visualize_significant_modules_as_networks(seurat_obj, DMEs, opt_data$results_path_ct_norm)
  	  # save DME csv - includes p-values, fold changes
  	  write.csv(DMEs, row.names=FALSE, quote=FALSE, file=paste0(opt_data$results_path_ct_norm, paste0("DMEs_",opt_test,".csv")))
    },
      error=function(cond){
      message("DME analysis failed!")
      message(cond)
      return(c())
    }
  )
}

