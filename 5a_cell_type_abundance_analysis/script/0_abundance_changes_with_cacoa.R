#author: Lisa Bast
#date: 13.09.2022
#about: cacoa for snRNAseq SCZ CTRL data (https://github.com/kharchenkolab/cacoa) to perform comparitive analysis of case vs. control
#       i) analysis of compositional shifts (cluster-based/ using CT annotation)
#       ii) analysis of transcriptional shifts (cluster-free/ not using CT annotation)
#       Fig. 1F
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
opt_cluster_based = TRUE #FALSE #
opt_run_per = "celltype" # "all" #"class" #
opt_inhibitory_only <- FALSE#TRUE #if TRUE runs for all inhibitory subtypes on the 37 cluster level
if (opt_cluster_based==FALSE){
  opt_graph_based <- TRUE
}
opt_species='human'
comp_name = 'longleaf'
opt_run_celltypes_in_parallel = TRUE
n_cores <- 35

code_path = getwd()
source("utils.R")

setwd("../")
project_path = getwd()
setwd("../")
main_path = getwd()

n_cluster=15 #37

if (opt_run_per=='celltype'){
  opt_show_cell_groups = FALSE
} else{
  opt_show_cell_groups = TRUE
}

if (opt_cluster_based==FALSE && opt_run_celltypes_in_parallel){
  args <- as.integer(commandArgs(trailingOnly=TRUE)) # ids referring to all cell types for running them in parallel
} else if (opt_run_per!='celltype' || opt_cluster_based==TRUE){
  args <- NaN
} else{
  args <- 13 #oligodendrocytes for testing
}

if (opt_inhibitory_only){
	results_path <- paste0(project_path,"results/Cell_type_abundance/Inhibitory_neurons/")
} else{
	results_path <- paste0(project_path,"results/Cell_type_abundance/all_cells/")
}
#make sure results path exists
if (file.exists(results_path)==FALSE) {
    dir.create(results_path, recursive=TRUE)
}

#Input Loom file(): Ct annotated and clustered cell type classes, downsampled to 20% of cells
data_path<-paste0(main_path,'/4_data_integration_and_cell_type_annotation/output/')

#get color_palettes:
palettes <- get_color_palettes(data_path)

#file all classes:
if (opt_inhibitory_only ==TRUE){
  filename <- "Samples_Inhibitory_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"
  n_cluster <- 37
} else{
  if (opt_run_per=='class'){
    CT_classes <- c("NonNeuronal","Excitatory","Inhibitory")
    filename <- paste0("Samples_",CT_classes[args],"_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom")
  } else{
    filename <- "Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"
    #ct_variable <- "cluster_name_15CTs"
  }
}


#load raw expression count matrix
data <- load_data_for_cacoa(data_path,filename)
print("loading data works")

if (opt_cluster_based){
  cao <- build_cacoa_object(data,n_cluster,n_cores,opt_cluster_based,args,palettes,opt_run_per)
  
  #close loom file connection
  data$close_all()
  
  cluster_based_analysis(cao, results_path)

} else{
  
  CT_id <- args[1] # is the input integer coming from .sh script in case of running on a server and the id specified above if this is just a test on a local machine
  #id refers to:
  #[1] "Excitatory Layer 2-3 IT neurons II"     "none (removed)"                         "Excitatory Layer 5-6 CT and NP neurons"
  #[4] "Excitatory Layer 3-4 IT neurons"        "Excitatory Layer 3-6 IT neurons"        "Excitatory Layer 2-3 IT neurons I"     
  #[7] "Excitatory Layer 5-6 IT neurons I"      "Excitatory Layer 5-6 IT neurons II"     "Inhibitory PVALB neurons"              
  #[10] "Inhibitory LAMP5 neurons"               "Inhibitory VIP neurons"                 "Inhibitory SST neurons"                
  #[13] "Oligodendrocytes"                       "Astrocytes"                             "Oligodendrocyte progenitor cells"      
  #[16] "Microglial cells"                       "Endothelial and mural cells" 
  cao <- build_cacoa_object(data,n_cluster,n_cores,opt_cluster_based,CT_id,palettes,opt_run_per)
  print("building cao object works")
  
  if (opt_run_per=='class'){
    CT_name <- CT_classes[CT_id]
  }else if (opt_run_per == 'all'){
    CT_name <- 'all'
  }else if (opt_run_per == 'celltype'){
    CT_name <- unique(data[[paste0(paste0('col_attrs/cluster_name_',n_cluster),'CTs')]][])[CT_id]
  }

  #close loom file connection
  data$close_all()
  #cluster-free compositional changes:
  #need to create Seurat object
  
  if (opt_graph_based){
    cao$estimateCellDensity(method='graph',n.cores=n_cores)
    print("cao$estimateCellDensity(method='graph',n.cores=n_cores) check!")
    
    cao$estimateDiffCellDensity(type='permutation',n.permutations=1000,verbose=TRUE,n.cores=n_cores)
    print(paste0("cao$estimateDiffCellDensity(type='wilcox',n.cores=n_cores) ",CT_name," check!"))
    
    gg_embedding <- cao$plotEmbedding(color.by='cell.groups')
    gg_z <- cao$plotDiffCellDensity(adjust.pvalues=FALSE,color.range=c("1%", "99%"),legend.position=c(0, 1))
    gg_z_adj <- cao$plotDiffCellDensity(adjust.pvalues=TRUE,color.range=c("1%", "99%"),legend.position=c(0, 1))
    
    plot_graph_based_results(cao,CT_name,results_path,gg_embedding,gg_z,gg_z_adj)
  } else{
    #plot density in SCZ and CTRL:
    cao$estimateCellDensity(method='kde',bins=400,name='cell.density.kde',n.cores=n_cores)
    print('cao$estimateCellDensity() successful')
    
    #cao$estimateDiffCellDensity(type='permutation', n.permutations=1000, verbose=TRUE, 
    #                            n.cores=n_cores, name='cell.density.kde')
    cao$estimateDiffCellDensity(type='permutation', n.permutations=200, verbose=TRUE, 
                                n.cores=n_cores, name='cell.density.kde')
    print('cao$estimateDiffCellDensity() successful')
    
    gg_embedding <- cao$plotEmbedding(color.by='cell.groups')
    gg_z_kde <- cao$plotDiffCellDensity(name='cell.density.kde',adjust.pvalues=FALSE,color.range=c("1%", "99%"),legend.position=c(0, 1),type='permutation')
    gg_z_adj_kde <- cao$plotDiffCellDensity(name='cell.density.kde',adjust.pvalues=TRUE,color.range=c("1%", "99%"),legend.position=c(0, 1),type='permutation')
    
    plot_density_based_results(cao,CT_name,results_path,opt_show_cell_groups,gg_embedding,gg_z_kde,gg_z_adj_kde)
  }
}

