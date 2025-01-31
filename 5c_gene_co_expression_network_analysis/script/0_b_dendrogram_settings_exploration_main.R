# about:  explore dendrogram settings
#         runs hdWGCNA (https://github.com/smorabit/hdWGCNA, https://www.biorxiv.org/content/10.1101/2022.09.22.509094v1.full.pdf),
#         some functions are slightly adapted to output additional values or change plotting functions, maths are not touched
# author: Hayley French, revised by Lisa Bast
# date: 4th Sept 2023
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

#define paths and load functions:
code_path = getwd()
source("utils_settings.R")
source("utils_hdWGCNA.R")
source("utils_plotting.R")

##### SETTINGS #####
num_threads = 24

#opt_save_obj <- FALSE
opt_run_gprofiler <- TRUE
opt_load_seurat_obj_with_metacells <- TRUE
opt_save_seurat_obj_with_metacells <- FALSE
opt_filter_genes <- "variable"#"fraction", "all", or "custom"

if (opt_load_seurat_obj_with_metacells == TRUE){
  library(foreach)
  # To Do: cell type id is input parameter
  CT_id <- commandArgs(trailingOnly=TRUE) # takes the input value specified in .sh file that runs this function
} else{
  CT_id <- 1
}

#get data settings
opt_data <- get_opt_data_hdWGCNA()
#get WGCNA settings
settings <- get_settings_hdWGCNA(opt_data)
settings$detectCutHeight = 0.975 # default is 0.995

#get paths
setwd("../")
setwd("../")
opt_data$main_path  = getwd()
opt_data <- get_paths_hdWGCNA(opt_data, settings, opt_filter_genes)

# set parameters to test
if (opt_data$bool_subsampled_data){
  deepSplit_range <- c(2,3,4) # simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive
  mergeCutHeight_range <- c(0.15, 0.2, 0.4) # dendrogram cut height for module merging
  minModuleSize_range <- c(25,50) # minimum module size for module detection
} else{
  deepSplit_range <- c(2,3,4) # simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive
  mergeCutHeight_range <- c(0.15, 0.2, 0.4) # dendrogram cut height for module merging
  minModuleSize_range <- c(25,50) # minimum module size for module detection
}


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

#create/ load seurat object with metacells
L <- get_metacells(opt_load_seurat_obj_with_metacells, opt_data, settings, opt_filter_genes)
seurat_obj <- L[[1]]
celltype_list <- L[[2]]
n_metacells_per_group <- L[[3]]
network_gene_list <- L[[4]]

i <- as.integer(CT_id)
ct <- as.character(celltype_list[i])
ct_ <- gsub(" ", "_", ct)
ct_ <- gsub("/", "-", ct_)

opt_data$results_path <- paste0(opt_data$results_path,"dendro_settings_expl/")
opt_data$results_path_ct_norm <- create_ct_subfolder(opt_data,ct_)
output_list <- list()
for (deepSplit in deepSplit_range){
  settings$deepSplit <- deepSplit
  for (mergeCutHeight in mergeCutHeight_range){
    settings$mergeCutHeight <- mergeCutHeight
    for (minModuleSize in minModuleSize_range) {
      settings$minModuleSize <- minModuleSize
      
      print(paste0("deepSplit: ",as.character(deepSplit)))
      print(paste0("mergeCutHeight: ",as.character(mergeCutHeight)))
      print(paste0("minModuleSize: ",as.character(minModuleSize)))
      
      exploration_settings_name <- paste0('deepSplit_',as.character(deepSplit), '_mergeCutHeight_', gsub('.', '', as.character(mergeCutHeight), fixed = TRUE),
                                           '_minModuleSize_', as.character(minModuleSize))
      #browser()
      if (opt_filter_genes=="custom"){
        #acvtivate a different set of genes but same metacells as before (created based on network genes for any of the cell types)
        seurat_obj <- SetupForWGCNA(
          seurat_obj,
          gene_select = opt_filter_genes,
          gene_list = network_gene_list[[ct]],
          wgcna_name = exploration_settings_name,
          metacell_location = "all_CTs_all_network_genes")
      } else{
        #acvtivate a different set of genes but same metacells as before (created based on network genes for any of the cell types)
        seurat_obj <- SetupForWGCNA(
          seurat_obj,
          gene_select = opt_filter_genes, #How to select genes? Select "variable", "fraction", "all", or "custom".
          wgcna_name = exploration_settings_name,
          metacell_location = "all_CTs_all_network_genes")
      }
      
      opt_data <- create_module_subfolder(opt_data,deepSplit,mergeCutHeight,minModuleSize)
      seurat_obj <- network_analysis_exploration(seurat_obj,ct,ct_,opt_data,settings,exploration_settings_name)
      if (length(seurat_obj)==0){
        next
      }
      # get active modules
      modules <- GetModules(seurat_obj, exploration_settings_name)
      modules = subset(modules, select = -c(module) )
      write.csv(modules, paste0(opt_data$results_path_long, 'modules.csv'))
      output_list[[exploration_settings_name]] <- modules
    }
  }
}

# merge the data frames in the list
merged_data <- bind_rows(output_list, .id = "Iteration")
merged_data <- merged_data %>%
  pivot_wider(names_from = "Iteration", values_from = "color")
write.csv(merged_data, paste0(opt_data$results_path_ct_norm, "modules_settings_comparison.csv"))
color_df <- merged_data %>% remove_rownames %>% column_to_rownames(var="gene_name")
#make sure inf and na values get removed (not sure if necessary, net was more the problem with the plot)
#keepers <- !sapply(color_df, function(i) any(sapply(i, is.infinite)))
#color_df <- color_df[keepers]
#keepers <- !sapply(color_df, function(i) any(sapply(i, is.na)))
#color_df <- color_df[keepers]

# get consensus dendrogram - use first combination as reference?
net <- GetNetworkData(seurat_obj)# do not adds "wgcna_name=exploration_settings_name", it will fail to load the network data and the graphic will produce an error

graphic_filename <- paste0(opt_data$results_path_ct_norm,"dendrogram_settings_comparison.pdf")
pdf(file = graphic_filename)#,   # The directory you want to save the file in
    #width = 12, # The width of the plot in inches
    #height = 8) # The height of the plot in inches
# plot dendrogram using WGCNA function
WGCNA::plotDendroAndColors(
  net$dendrograms[[1]],
  color_df,
  groupLabels=colnames(color_df),
  dendroLabels = FALSE,
  hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "Module parameters exploration dendrogram",
  marAll = c(0.5, 16, 2, 0.05) # bottom, left, top, right
)

dev.off()

