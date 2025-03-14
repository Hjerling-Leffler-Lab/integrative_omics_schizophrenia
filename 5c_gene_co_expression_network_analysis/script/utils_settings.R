## about: settings for running hdWGCNA
## author: Lisa Bast
## date: 29.06.2023
## version: 0.0.1

#get data settings
get_opt_data_hdWGCNA <- function(){
  opt_data <- c()
  
  #which data should be extracted from seurat object:
  opt_data$bool_subsampled_data = FALSE

  # define cell type level
  opt_data$level <- "cluster_name_15CTs"
  
  #for seurat object:
  opt_data$feature_name = "Gene"
  opt_data$cell_name = "CellID"
  
  opt_data$sample_variable <- "Donor"
  opt_data$group <- "Disease"

  opt_data$normalize <- "basic" #"counts" #"metacells" # normalization with sctransform
  #basic: runs NormalizeData() and later also sclaeData() and FindVariableFeatures() on raw counts and performs basic noralization of metacells
  #counts: runs sctransform() on raw counts while regressing certain varibles out and performs basic noralization of metacells
  #metacells: runs NormalizeData() on raw counts and sctransform() on metacells while regressing certain varibles out
  opt_data$bool_load_normalized_data = FALSE # if TRUE and opt_data$normalize=="counts" it directly loads the already sctransform normalized counts
  
  if (opt_data$normalize=="basic"){
    opt_data$vars_to_regress <- NULL
  } else{
    opt_data$vars_to_regress <- c("Donor")
  }
  
  return(opt_data)
}

#get settings for algorithm:
get_settings_hdWGCNA <- function(opt_data){
  settings <- c()
  
  settings$randomSeed = 8747
  
  settings$gene_sign_score_comp <- 'Seurat'#'UCell' # UCell package requires R >= 4.2.0, make sure it is installed

  settings$export_TOM_matrix <- TRUE # takes pretty long
  settings$fraction_of_cells_expr_gene <- 0.05 # initial data filtering before metacells are generated
  settings$dim_reduction='pca'#'harmony'
  settings$wgcna_name = 'scz'
  # traits for module trait correlation analysis
  settings$traits <- c('Disease', 'Age', 'Sex', 'PMI_h', 'nCount_RNA','nFeature_RNA','pct_counts_mito','pct_counts_ribo')
  
  if (opt_data$bool_subsampled_data==TRUE){
    settings$meta_k <- 4#Number of nearest neighbors to aggregate. Default = 50
    settings$meta_max_shared <- 1#the maximum number of cells to be shared across two metacells
    settings$meta_min_cells <- 16# the minimum number of cells in a particular grouping to construct metacells
    n_cells_dataset <- 43000
  } else{
    settings$meta_k <- 30 #50 old setting #Number of nearest neighbors to aggregate. Default = 50
    settings$meta_max_shared <- 6#13 old setting #the maximum number of cells to be shared across two metacells
    settings$meta_min_cells <- 120#250 old setting #400 was too strict, too many cell types removed# the minimum number of cells in a particular grouping to construct metacells
    n_cells_dataset <- 430000
  }
  settings$target_number_meta_cells <- round(n_cells_dataset/settings$meta_k)
  
  if (opt_data$normalize == "counts"){
    settings$slot_sel <- "scale.data"
    settings$assay_sel <- "SCT"
    settings$wgcna_name <- paste0(settings$wgcna_name,"_SCT_",opt_data$normalize)
  } else if (opt_data$normalize == "metacells"){
    settings$slot_sel <- "counts"
    settings$assay_sel <- "RNA"
    settings$wgcna_name <- paste0(settings$wgcna_name,"_SCT_",opt_data$normalize)
  } else{ # "basic" --> keep defaults
    settings$slot_sel <- "counts"
    settings$assay_sel <- "RNA"
  }

  ### universal WGCNA settings:
  # module detection settings
  settings$opt_use_power_12 <- FALSE
  if (settings$opt_use_power_12) {
    settings$power_str <- "fixed_power_12"
  } else {
    settings$power_str <- "power_detected"
  }
  settings$corType = "bicor" #"pearson" # for gene-to-gene correlation, but what implication has bicor (bidweigth midcorrelation) actually?
  settings$TOMType = "signed"
  settings$networkType = "signed"
  
  #these settings might be updated later in update_dendrogram_settings() function !!!!
  settings$mergeCutHeight = 0.2 #new: 0.1 # default is 0.2, in pseudobulk, rather high value of mergeCutHeight seems to be beneficial; for fewer samples, higher mergeCutHeight; for larger modules, higher mergeCutHeight
  settings$detectCutHeight = 0.995 # default is 0.995
  settings$deepSplit = 4 # For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. The higher the value (or if TRUE), the more and smaller clusters will be produced. Default for hdWGCNA is 4.
  settings$minModuleSize = 25 # default is 50, minimum number of genes per module

  # define if MEs should be scaled (MEs = module eigengenes)
  
  return(settings)
}

#define data_counts_filename
get_data_counts_filename <- function(opt_data){
  if (opt_data$bool_load_normalized_data){
    norm_str <- "_sct_norm_v1"
  } else{
    norm_str <- ""
  }
    
  if (opt_data$bool_subsampled_data==TRUE){
    subsample_str <- "_subsampled_to_10_percent_cells"
  } else{
    subsample_str <- ""
  }
  data_counts_filename <- paste0("Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered",subsample_str,norm_str,".loom")

  
  return(data_counts_filename)
}



#get paths
# define paths:
get_paths_hdWGCNA <- function(opt_data, settings, opt_filter_genes){

  .libPaths( c( .libPaths(), "") )
  
  opt_data$code_path <- paste0(opt_data$main_path,"/code/5c_gene_co_expression_network_analysis/")
    
  opt_data$data_counts_path <- paste0(opt_data$main_path,"/code/4_data_integration_and_cell_type_annotation/output/")
  opt_data$data_background_path <- paste0(opt_data$code_path,"data/all_genes_detected_in_snRNAseq.csv")
  if (opt_data$bool_subsampled_data){
    add_folder = "subsampled_data/"
  }else{
    add_folder = ""
  }
  if (settings$wgcna_name == 'scz_metacell_test'){
    opt_data$results_path <- paste0(opt_data$code_path,"output/", opt_data$level, "/", settings$wgcna_name, "/",add_folder)
  } else{
    opt_data$results_path <- paste0(opt_data$code_path,"output/", opt_data$level, "/", settings$wgcna_name, "/cor_", settings$corType,"/",add_folder)
  }

  opt_data$data_counts_filename <- get_data_counts_filename(opt_data)
  opt_data$results_path <- paste0(opt_data$results_path,opt_filter_genes,"/")

  #make sure results path exists
  if (file.exists(opt_data$results_path)==FALSE) {
    dir.create(opt_data$results_path, recursive=TRUE)
  }
  return(opt_data)
}

get_settings_metacell_exploration <- function(opt_data){
  
  settings <- c()
  settings$fraction_of_cells_expr_gene <- 0.05
  settings$traits <- c('Disease', 'Age', 'Sex', 'PMI_h', 'nCount_RNA','nFeature_RNA','pct_counts_mito','pct_counts_ribo')
  settings$randomSeed = 8747
  settings$export_TOM_matrix <- FALSE # takes pretty long
  settings$gene_sign_score_comp <- 'Seurat'#'UCell' # UCell package requires R >= 4.2.0, make sure it is installed
  settings$dim_reduction='pca'#'harmony'
  settings$wgcna_name = 'scz_metacell_test'
  if (opt_data$normalize == "counts"){
    settings$slot_sel <- "scale.data"
    settings$assay_sel <- "SCT"
  } else if (opt_data$normalize == "metacells"){
    settings$slot_sel <- "counts"
    settings$assay_sel <- "RNA"
  } else{ # "basic" --> keep defaults
    settings$slot_sel <- "counts"
    settings$assay_sel <- "RNA"
  }
  
  return(settings)
}

update_dendrogram_settings <- function(settings, ct, opt_filter_genes){
  if (opt_filter_genes=='variable'){
    settings$minModuleSize = 25 
    settings$deepSplit = 2 #default for most CTs
    #all cell type names:
    #"Excitatory Layer 3-6 IT neurons"        "Excitatory Layer 2-3 IT neurons II"     "Excitatory Layer 2-3 IT neurons I"     
    #"Excitatory Layer 3-4 IT neurons"        "Excitatory Layer 5-6 CT and NP neurons" "Excitatory Layer 5-6 IT neurons I"     
    #"Excitatory Layer 5-6 IT neurons II"     "Inhibitory PVALB neurons"               "Inhibitory LAMP5 neurons"              
    #"Inhibitory VIP neurons"                 "Inhibitory SST neurons"                 "Oligodendrocyte progenitor cells"      
    #"Microglial cells"                       "Oligodendrocytes"                       "Astrocytes"                            
    #"Endothelial and mural cells"  
    if (ct %in% list("Excitatory Layer 2-3 IT neurons I" , "Inhibitory PVALB neurons" , "Oligodendrocytes", "Oligodendrocyte progenitor cells" ,"Excitatory Layer 5-6 CT and NP neurons")){
      settings$deepSplit = 3
    }
    settings$mergeCutHeight = 0.15 #default for most CTs
    if (ct %in% list("Astrocytes","Excitatory Layer 5-6 IT neurons I" ,"Inhibitory PVALB neurons" )){
      settings$mergeCutHeight = 0.2 
    } else if (ct %in% list("Excitatory Layer 2-3 IT neurons I","Excitatory Layer 2-3 IT neurons II")){
      settings$mergeCutHeight = 0.4
    }
  }
  return(settings)
}