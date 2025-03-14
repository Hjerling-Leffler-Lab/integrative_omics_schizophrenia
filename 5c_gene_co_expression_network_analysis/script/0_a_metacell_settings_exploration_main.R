# about: explore metacell settings, run for different hyperparameters to obtain results for all opt_filter_genes options
#        runs hdWGCNA (https://github.com/smorabit/hdWGCNA, https://www.biorxiv.org/content/10.1101/2022.09.22.509094v1.full.pdf),
#        some functions are slightly adapted to output additional values or change plotting functions, maths are not touched
# author: Lisa Bast, revised by Hayley French
# date: 15th June 2023
# version: 0.0.1

library(Seurat)
library(SeuratDisk)
library(harmony)
library(tidyverse)
library(cowplot)
library(patchwork)
library(hdWGCNA)
library(WGCNA)
library(plyr)
library(ggplot2)
library(pracma)#for ceil function
library(stringr) # for str_replace()
library(foreach)

#define paths and load functions:
code_path = getwd()
source("utils_settings.R")
source("utils_hdWGCNA.R")
source("utils_plotting.R")


##### General SETTINGS #####
num_threads = 24

calculate_full_matrix_sparsity <- FALSE
opt_calculate_correlation <- FALSE # computationally intense (expensive on memory)
opt_load_network_gene_list <- FALSE


# cell type id is input parameter
opt_filter_genes_id <- as.integer(commandArgs(trailingOnly=TRUE)) # takes the input value specified in .sh file that runs this function

opt_filter_genes_list <- c("fraction","variable", "all", "custom")
opt_filter_genes <- opt_filter_genes_list[opt_filter_genes_id]

#get data settings
opt_data <- get_opt_data_hdWGCNA()

settings <- get_settings_metacell_exploration(opt_data)
if (opt_filter_genes=="custom"){
  add_str <- paste0("_fraction_genes_", str_replace(as.character(settings$fraction_of_cells_expr_gene),"[.]",""))
} else{
  add_str <- ""
}

if (opt_data$bool_subsampled_data){
  meta_min_cells_range <- c(8,12,16)# the minimum number of cells in a particular grouping (per cell) to construct metacells
  meta_k_range <- c(4,6,8) #Number of nearest neighbors to aggregate. Default = 50
  meta_max_shared_range <- c(1,2,3) # should have the same length as meta_min_cells_range
} else{
  meta_min_cells_range <- c(80,100,120,250) # the minimum number of cells in a particular grouping (per cell) to construct metacells
  meta_k_range <- c(30,50,70) #Number of nearest neighbors to aggregate. Default = 50
  meta_max_shared_range <- ceil(0.05*meta_min_cells_range) #20#the maximum number of cells to be shared across two metacells
  # meta_max_shared_range should have the same length as meta_min_cells_range
}

#get paths
setwd("../")
setwd("../")
opt_data$main_path  = getwd()
opt_data <- get_paths_hdWGCNA(opt_data, settings, opt_filter_genes)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
if (num_threads>1){
  enableWGCNAThreads(nThreads = num_threads)
}

seurat_obj <- get_seurat_object_ready(settings,opt_data)
n_cells_dataset <- dim(seurat_obj)[2]
#some preparation before metacells can be generated:
seurat_obj <- metacell_preparations(seurat_obj,settings,opt_data)

DF_R <- visualize_number_cells_per_donor_proportion(seurat_obj,opt_data, add_str)

if (opt_calculate_correlation){
  colors <- c("orange","blue","green","violet","black","darkgrey")
  counter <- 1
}

if (opt_filter_genes=="custom"){
  
  if (opt_load_network_gene_list==TRUE){
    network_gene_list <- readRDS(paste0(opt_data$results_path, "network_gene_list.RData"))
  } else{ 
    # get list of genes for each cell type
    network_gene_list <- select_genes_for_network_analysis(seurat_obj, fraction = settings$fraction_of_cells_expr_gene, group.by = opt_data$level)
    saveRDS(network_gene_list, file=paste0(opt_data$results_path, "network_gene_list.RData"))
  }
  all_network_genes <- unique(unlist(network_gene_list))
  display_genes_filtering_stats(celltype_list,all_network_genes,network_gene_list,seurat_obj,settings$fraction_of_cells_expr_gene)
} 

# loop over arrays of hyper parameter values:
for (k in meta_k_range){
  settings$meta_k <- k
  settings$target_number_meta_cells <- round(n_cells_dataset/settings$meta_k)
  for (m in meta_min_cells_range){
    settings$meta_min_cells <- m
    i <- which(meta_min_cells_range==m)
    settings$meta_max_shared <- meta_max_shared_range[i]
    print(paste0("Number neighbors: ",as.character(k)))
    print(paste0("Min cells: ",as.character(m)))
    print(paste0("Max cells shared: ",as.character(settings$meta_max_shared)))
    if (opt_calculate_correlation){
      color <- colors[counter]
    }
    opt_data <- create_subfolder(opt_data,k,settings$meta_max_shared,m)
    
    #sparsity for full matrix:
    if ((k==meta_k_range[1]) & (m==meta_min_cells_range[1])){
      #initialize sparsity matrix:
      n_donors_ctrl <- length(unique(seurat_obj$Donor[seurat_obj$Disease=='CTRL']))
      n_donors_case <- length(unique(seurat_obj$Donor[seurat_obj$Disease=='SCZ']))
      n_CTs <- length(unique(seurat_obj$cluster_name_15CTs))
      
      M_full <- GetAssayData(seurat_obj,'counts')
      if (calculate_full_matrix_sparsity == TRUE) {
        S_R <- data.frame(k_nns = 0, cells_max_shared = Inf, cells_min = 0, sparsity=calculate_sparsity(M_full), number_of_donors_ctrl=n_donors_ctrl, number_of_donors_case=n_donors_case, number_cell_types = n_CTs)
      } else {
        S_R <- data.frame(k_nns = 0, cells_max_shared = Inf, cells_min = 0, sparsity=NA, number_of_donors_ctrl=n_donors_ctrl, number_of_donors_case=n_donors_case, number_cell_types = n_CTs)
      }
    }
    
    #generate metacells by grouping by Sample and cell_type to achieve the desired result.
    #exclude underrepresented cell types with min_cells
    if (opt_filter_genes=="custom"){
      #determine which slots in seurat_obj are used:
      seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = opt_filter_genes,
        gene_list = all_network_genes,#network_gene_list[[ct]],#shouldn't we have genes here expressed in at least one cell type and after metacell construction run SetupForWGCNA() for this Cell types genes?
        wgcna_name = settings$wgcna_name)
    } else{
      # only run metacell aggregation one time
      print('Setting up and running metacells.')
      
      #determine which slots in seurat_obj are used:
      seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = opt_filter_genes, # upcoming: explore other settings further
        wgcna_name = settings$wgcna_name)
    }
    
    
    #generate metacells by grouping by Sample and cell_type to achieve the desired result.
    #exclude underrepresented cell types with min_cells
    seurat_obj_updated <- tryCatch(
      {
        seurat_obj_updated <- MetacellsByGroups(
          seurat_obj = seurat_obj,
          group.by = c(opt_data$level, opt_data$sample_variable, opt_data$group), # specify the columns in seurat_obj@meta.data to group by
          ident.group = opt_data$level,  # set the Idents of the metacell seurat object
          reduction = settings$dim_reduction, # select the dimensionality reduction to perform KNN on
          k = settings$meta_k, # nearest-neighbors parameter
          max_shared = settings$meta_max_shared, # maximum number of shared cells between two metacells
          min_cells = settings$meta_min_cells,
          mode='average',
          slot = settings$slot_sel,
          assay = settings$assay_sel,
          target_metacells = settings$target_number_meta_cells,
          wgcna_name=settings$wgcna_name
        )
      },
      error=function(cond){
        message("MetacellsByGroups() failed!")
        message(cond)
        return(NA)
      }
    )
    if (is.na(seurat_obj_updated)){
      n_donors_ctrl <- 0
      n_donors_case <- 0
      n_CTs <- 0
      s_r <- data.frame(k_nns = k, cells_max_shared = settings$meta_max_shared, cells_min = m, sparsity=NA, number_of_donors_ctrl=n_donors_ctrl, number_of_donors_case=n_donors_case, number_cell_types = n_CTs)
      S_R <- rbind(S_R,s_r)
      next
    } else{
      seurat_obj <- seurat_obj_updated
    }
    
    #calculate sparsity of metacell matrix
    metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name = settings$wgcna_name)
    M <- GetAssayData(object = metacell_obj, slot = "counts")

    #which and how many donors are contributing to metacell:
    list_groups_contributing_to_metacells <- unique(vapply(strsplit(colnames(M),'_',fixed=TRUE),  "[","",1))
    #dput(list_groups_contributing_to_metacells, paste0(opt_data$results_path_long,"list_groups_contributing_to_metacells.txt"))
	  dput(list_groups_contributing_to_metacells, paste0(opt_data$results_path_long,"list_groups_contributing_to_metacells_",add_str,".txt"))

    #compute, store and plot sparsity for count matrix and metacell matrices for a range of parameters
    n_donors_ctrl <- length(unique(metacell_obj$Donor[metacell_obj$Disease=='CTRL']))
    n_donors_case <- length(unique(metacell_obj$Donor[metacell_obj$Disease=='SCZ']))
    n_CTs <- length(unique(metacell_obj$cluster_name_15CTs))

    s_r <- data.frame(k_nns = k, cells_max_shared = settings$meta_max_shared, cells_min = m, sparsity=calculate_sparsity(M), number_of_donors_ctrl=n_donors_ctrl, number_of_donors_case=n_donors_case, number_cell_types = n_CTs)
    S_R <- rbind(S_R,s_r)
    #save sparsity result
    write.csv(S_R, paste0(opt_data$results_path,"sparsity_result",add_str,".csv"))
    
    if (opt_calculate_correlation){
      #plot pairwise gene correlation distribution for range of settings
      if ((k==meta_k_range[1]) && (m==meta_min_cells_range[1])){
        corr_vec_full <- as.vector(sparse.cor3(Transpose(M_full)))
        d_full <- density(corr_vec_full)
        density_plot <- densityPlot(d_full, histo="none", xlim=c(-1,1), lwd = 2, col = "red", main = "Density of row correlations", label = "full count matrix")
      }
      corr_vec <- as.vector(sparse.cor3(Transpose(M)))
      density_plot <- density_plot + densityplot(corr_vec, histo="none", xlim=c(-1,1), col=color, lwd=2, label = paste0(as.character(k),'_neighbors_',as.character(settings$meta_max_shared),'_cells_shared_at_least_',as.character(m),'_cells'))
      if ((k==meta_k_range[length(meta_k_range)]) && (m==meta_min_cells_range[length(meta_min_cells_range)])){
        graphic_filename <- paste0(opt_data$results_path,"row_correlation_desities.pdf")
        pdf(file = graphic_filename,   # The directory you want to save the file in
            width = 12, # The width of the plot in inches
            height = 9) # The height of the plot in inches
        print(density_plot)
        dev.off()
      }
      counter <- counter+1
    }
    
    #normalize metacell expression matrix:
    #normalize metacell expression matrix:
    if (opt_data$normalize == "basic"){
      seurat_obj <- NormalizeMetacells(seurat_obj)
    } else if (opt_data$normalize == "metacells"){
      metacell_obj <- GetMetacellObject(seurat_obj)
      meta_data <- seurat_obj@meta.data[,-1] # gets everything that is stored under metadata except orig.ident
      metacell_obj <- normalize_data_with_sctransform(metacell_obj, opt_data)
      # only keep genes that were used for SCTransform
      sct_data <- GetAssayData(metacell_obj, slot='scale.data', assay='SCT')
      sct_genes <- rownames(sct_data)
      gene_list <- GetWGCNAGenes(seurat_obj)
      gene_list <- gene_list[gene_list %in% sct_genes]
      # update the genes used for WGCNA, and reset the metacell object
      seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list)
      seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj)
    }

    #visualize metacells:
    #needs trouble shooting
    tryCatch(
      visualize_meta_cells(seurat_obj, settings, opt_data),
      error=function(cond){
        message("visualize_meta_cells() failed!")
        message(cond)
        return(NA)
      }
    )
    
  }
}


