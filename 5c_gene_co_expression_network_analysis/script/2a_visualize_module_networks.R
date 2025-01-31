# about: network visualization for specific cell types and modules
#        for each cell type plot only significant or only bordeaux modules
#        highlight genes based on different criteria (mBATcombo, DAPs, ...)
#        Fig. 5C
# author: Lisa Bast 
# date:  28th Februari 2024
# version: 0.0.1

library("ggraph")
library("igraph")
library("tidygraph")
library("Seurat")
library("SeuratDisk")

code_path = getwd()
source("utils_plotting.R")

##### SETTINGS #####
opt_load_from_seurat_obj <- TRUE#FALSE
setwd("../")
setwd("../")
main_path  = getwd()

opt_data <- c()

#which data should be extracted from seurat object:
opt_data$bool_subsampled_data = FALSE
code_path <- paste0(main_path,"/code/5c_gene_co_expression_network_analysis/")
results_path <- paste0(code_path,"output/cluster_name_15CTs/scz/cor_bicor/")
TOM_path <- paste0(code_path,"TOM/") 
seurat_obj_path <- paste0(code_path,"seurat_objects/")
GM_path <- paste0(main_path,"/5b_differential_gene_expression_analysis/data/")
metadata_path <- paste0(main_path,"data/")
                                  
highlight_colors <- c("turquoise","purple","blue","orange","red")
other_color <- "black"

## networks with genes highlighted
#a) trubeskoy genes from gene matrix (scz2022PriorityGene)
GM <- read.csv(paste0(GM_path,"geneMatrix_v2.tsv"),sep = "\t",header = TRUE)
WES_genes <- GM$gene_name[GM$scz2022.wes.qvalSig05==TRUE]
common_variants <- GM$gene_name[GM$scz2022PriorityGene==TRUE]
#b) mbatcombo genes
mBAT_genes <- read.csv(paste0(metadata_path,"T14_mBATcombo_bordeaux_module_genes.tsv"),sep = "\t",header = TRUE)
colnames(mBAT_genes)[colnames(mBAT_genes) == "ENSGID"] ="ensgid"
gene_ensid_mapping <- read.csv(paste0(GM_path,"Gene_ensgid_mapping.csv"),header = TRUE)
mBAT_genes <- merge(mBAT_genes,gene_ensid_mapping,by="ensgid",all.x=TRUE)
mBAT_genes <- mBAT_genes[mBAT_genes$FDR_mBATcombo_onlyBordeauxGenes<0.05,]
#c) DAP genes
DAPs <- read.csv(paste0(metadata_path,"sign_DAP_genes.csv"),header=FALSE)
DAP_genes <- DAPs$V1

#files <- list.files(TOM_path)
files <- c("Excitatory_Layer_5-6_IT_neurons_I__norm_basic_TOM.Rds",               
           "Inhibitory_LAMP5_neurons__norm_basic_TOM.Rds",                         
           "Inhibitory_PVALB_neurons__norm_basic_TOM.Rds",                         
           "Inhibitory_VIP_neurons__norm_basic_TOM.Rds",                          
           "Oligodendrocytes__norm_basic_TOM.Rds",
           "Inhibitory_SST_neurons__norm_basic_TOM.Rds",   
           "Astrocytes__norm_basic_TOM.Rds",                                        
           "Excitatory_Layer_2-3_IT_neurons_I__norm_basic_TOM.Rds",                
           "Excitatory_Layer_2-3_IT_neurons_II__norm_basic_TOM.Rds",               
           "Excitatory_Layer_3-6_IT_neurons__norm_basic_TOM.Rds",                   
           "Excitatory_Layer_5-6_CT_and_NP_neurons__norm_basic_TOM.Rds")          
           
for (file in files){
  print(file)
  if (endsWith(file,"_TOM.Rds")){
    ct = strsplit(file,"__")[[1]][1]
    print(ct)
    results_path_ct <- paste0(results_path,ct,"/norm_basic/")
    if (opt_load_from_seurat_obj){
      # load seurat_obj
      browser()
      seurat_obj <- LoadH5Seurat(paste0(seurat_obj_path,"S_",ct,".h5Seurat"))
      # get modules and TOM from the seurat obj
      modules <- GetModules(seurat_obj) %>% 
        subset(module != 'grey') %>% 
        mutate(module = droplevels(module))
      mods <- levels(modules$module)
      TOM <- GetTOM(seurat_obj)
    } else{
      TOM_file <- paste0(ct,"__norm_basic_TOM.Rds")
      module_file <- paste0(ct,"__norm_basic_modules_for_TOM.Rds")
      TOM <- readRDS(paste0(TOM_path,TOM_file))
      modules <- readRDS(paste0(TOM_path,module_file))
    }
    DMEs <- read.csv(paste0(results_path_ct,"DMEs_LR.csv"))
    DMEs_filt <- DMEs[DMEs$p_val_adj<=0.05,]
    
    if (opt_load_from_seurat_obj){
      visualize_significant_modules_as_networks(seurat_obj, DMEs, results_path_ct)
    }
    #filter only for one sign module at a time
    for (sign_module in DMEs_filt$module){
      module_sel <- modules[modules$module==sign_module,]
      TOM_sel <- TOM[module_sel$gene_name,module_sel$gene_name]
      mBAT_genes_sel <- mBAT_genes$gene_name[mBAT_genes$Celltype_Module==paste0(ct,"_",sign_module)]
      gene_lists <- list(common_variants = common_variants, WES_genes = WES_genes, mBAT_genes = mBAT_genes_sel, DAP_genes=DAP_genes)
      for (gene_list_name in c("all","common_variants","WES_genes","mBAT_genes","DAP_genes")){
        plot_module_network_with_genes_highlighted(TOM_sel,module_sel,gene_lists,gene_list_name,ct,sign_module,results_path_ct,highlight_colors,other_color)
      }
    }
  }
}



