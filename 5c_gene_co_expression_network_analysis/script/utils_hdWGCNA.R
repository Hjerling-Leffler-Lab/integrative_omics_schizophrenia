### Functions to load and set up data as seurat object and compute metacells
get_seurat_object_ready <- function(settings,opt_data, list_ensgid_of_interest = c()){
  # load seurat object
  filename_w_dir <- paste0(opt_data$data_counts_path,opt_data$data_counts_filename)
  data <- Connect(filename_w_dir)
    #sc <- import("scanpy")
    #adata <- sc$read_h5ad(filename_w_dir)
    #srt <- adata_to_srt(adata)
  seurat_obj <- as.Seurat(data, features = opt_data$feature_name, cells = opt_data$cell_name)
  
  #subset seurat if necessary
  if (length(list_ensgid_of_interest)>0){
    #first filter genes
    genes<-data[['row_attrs/Gene']][]
    gene_IDs <- data[['row_attrs/Accession']][]
    ID_sel <- which(gene_IDs %in% list_ensgid_of_interest) #[genes[i] for i in gene_IDs if i in list_ensgid_of_interest]
    genes_sel <- genes[ID_sel]
    seurat_obj <- subset(seurat_obj,features=genes_sel)
  }
  
  if (opt_data$opt_exclude_MT_genes){
    #get features of seurat obj
    genes <- rownames(seurat_obj)
    #remove MT genes
    ID_sel <- which(!startsWith(genes,"MT-"))
    seurat_obj <- subset(seurat_obj,features=genes[ID_sel])
  } 
  rm(data)
  
  #remove the "none (removed)" cluster
  seurat_obj <- subset(x=seurat_obj, subset=cluster_name_15CTs!="none (removed)")
  # Calculate the proportion of transcripts mapping to mitochondrial genes
  #add info about 'pct_counts_mt'
  seurat_obj[["pct_counts_mito"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj[["pct_counts_ribo"]] <- PercentageFeatureSet(object = seurat_obj, pattern="^RP[LS]")
	#make trait vars numerical
	if (opt_data$scdrs_variable %in% settings$traits){
		num_vars <- c('Age','PMI_h','nCount_RNA','nFeature_RNA','pct_counts_mito','pct_counts_ribo',opt_data$scdrs_variable)
	} else{
		num_vars <- c('Age','PMI_h','nCount_RNA','nFeature_RNA','pct_counts_mito','pct_counts_ribo')
	}
	seurat_obj <- make_seurat_var_numerical(num_vars, seurat_obj)
	#make trait vars factors 
	cat_vars <- c('Disease', 'Sex')
	seurat_obj <- make_seurat_var_factor(cat_vars, seurat_obj)
  
  #optionally normalize with sctransform
  if (opt_data$normalize == "counts"){
    meta_data <- seurat_obj@meta.data[,-1] # gets everything that is stored under metadata except orig.ident
    seurat_obj <- normalize_data_with_sctransform(seurat_obj, opt_data)
    feature_sel <- VariableFeatures(seurat_obj) # only the genes available in scale.data slot
  } else if (opt_data$normalize == "basic"){
    feature_sel <- NULL # keep default
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  } else if (opt_data$normalize == "metacells"){
    feature_sel <- NULL
  }
  
  #basic filtering: select genes that are expressed in at least 5% of cells in this dataset, and we will name our hdWGCNA experiment "scz".
  #seurat_obj <- SetupForWGCNA(
    #seurat_obj,
    #features = feature_sel,
    #gene_select = "fraction", # the gene selection approach
    #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    #wgcna_name = ct # the name of the hdWGCNA experiment
  #)
  
  return(seurat_obj)
}

make_seurat_var_numerical <- function(num_vars, seurat_obj){
  for (i in seq(1,length(num_vars))){
    v <- num_vars[i]
	  eval(parse(text=paste0("seurat_obj$",v," <- as.numeric(seurat_obj$",v,")")))
	}
	return(seurat_obj)
}

make_seurat_var_factor <- function(cat_vars, seurat_obj){
  for (i in seq(1,length(cat_vars))){
    v <- cat_vars[i]
	  eval(parse(text=paste0("seurat_obj$",v," <- as.factor(seurat_obj$",v,")")))
	}
	return(seurat_obj)
}

#normalize single cell data with sctransform, optionally called from get_seurat_object_ready
normalize_data_with_sctransform <- function(data, opt_data){
  if (opt_data$normalize=="counts"){
    data_norm <- SCTransform(data, 
                           verbose = FALSE,
                           latent_var = c("log_umi"), 
                           vars.to.regress = opt_data$vars_to_regress, 
                           return_gene_attr = TRUE) 
  } else if (opt_data$normalize=="metacells"){
    data_norm <- SCTransform(data, 
                            verbose = FALSE,
                            latent_var = c("log_umi"), 
                            vars.to.regress = opt_data$vars_to_regress, 
                            return.only.var.genes = FALSE) 
  }
  
  
  return(data_norm)
}

metacell_preparations <- function(seurat_obj,settings,opt_data){
  if (settings$dim_reduction == 'pca'){ #option a: PCA
    if (opt_data$normalize == "counts"){ # determine which genes to run PCA on
      feature_sel <- VariableFeatures(object = seurat_obj, 
                                      assay = settings$assay_sel)
    } else{
      #To DO: add option to regress out NormalizeData
      #double check: does this use slot data, assay RNA?
      seurat_obj <- ScaleData(seurat_obj,vars.to.regress = opt_data$vars_to_regress) # if opt_data$vars_to_regress !=NULL --> regressed out
      seurat_obj <- FindVariableFeatures(seurat_obj, 
                                         selection.method = 'vst', 
                                         assay = settings$assay_sel)
      #feature_sel <- VariableFeatures(object = seurat_obj, 
      #                                selection.method ='vst',
      #                                assay = settings$assay_sel)
      feature_sel <- VariableFeatures(object = seurat_obj, 
                                      method ='vst',
                                      assay = settings$assay_sel)
    }
    seurat_obj <- RunPCA(seurat_obj, 
                         npcs = 40, 
                         features = feature_sel, 
                         assay.type=settings$assay_sel)

  } else if (settings$dim_reduction == 'harmony'){ #option b: harmony data integration 
    seurat_obj <- RunHarmony(seurat_obj,
                             group.by.vars = opt_data$sample_variable, 
                             plot_convergence = TRUE, 
                             assay.use = settings$assay_sel)
  }
  seurat_obj <- RunUMAP(seurat_obj, 
                        reduction=settings$dim_reduction, 
                        dims=1:15, 
                        assay = settings$assay_sel, 
                        slot = settings$slot_sel)
  
  return(seurat_obj)
}

define_cell_type_markers <- function(seurat_obj,opt_data){
  if (file.exists(paste0(opt_data$results_path, "cell_type_markers.csv"))){
    # read in output file from Seurat FindAllMarkers()
    markers <- read.csv(paste0(opt_data$results_path, "cell_type_markers.csv"))
  } else { #compute cell type marker genes for module overlap
    Idents(seurat_obj) <- seurat_obj@meta.data[[opt_data$level]]
    markers <- Seurat::FindAllMarkers(
      seurat_obj,
      only.pos = TRUE,
      logfc.threshold=1
    )
    # save cell type markers
    write.csv(markers, row.names=FALSE, quote=FALSE, file=paste0(opt_data$results_path, "cell_type_markers.csv"))
  }
  return(markers)
}

select_genes_for_network_analysis <- function (seurat_obj, fraction, group.by) {
  assay <- DefaultAssay(seurat_obj)
  if (packageDescription("Seurat")$Version >= package_version("5.0")){
    expr_mat <- LayerData(seurat_obj, slot = "counts")
  } else{
    expr_mat <- GetAssayData(seurat_obj, slot = "counts")
  }
  n_chunks <- ceiling(ncol(expr_mat)/10000)
  if (n_chunks == 1) {
    chunks <- factor(rep(1), levels = 1)
  }
  else {
    chunks <- cut(1:nrow(expr_mat), n_chunks)
  }
  expr_mat <- do.call(rbind, lapply(levels(chunks), function(x) {
    cur <- expr_mat[chunks == x, ]
    cur[cur > 0] <- 1
    cur
  }))
  group_gene_list <- list()
  if (!is.null(group.by)) {
    groups <- unique(seurat_obj@meta.data[[group.by]])
    for (cur_group in groups) {
      cur_expr <- expr_mat[, seurat_obj@meta.data[[group.by]] == 
                             cur_group]
      print(dim(cur_expr))
      gene_filter <- Matrix::rowSums(cur_expr) >= 
        round(fraction * ncol(cur_expr))
      group_gene_list[[cur_group]] <- rownames(seurat_obj)[gene_filter]
    }
  } 
  return(group_gene_list)
}

get_metacells <- function(opt_load_seurat_obj_with_metacells, opt_data, settings, opt_filter_genes){
  network_gene_list <- c()
  #create/ load seurat object with metacells
  if (opt_load_seurat_obj_with_metacells==TRUE){
    tic("loading metacells total")
    #load seurat object with metacells
    tic("loading s_obj_w_metcells.rds")
    seurat_obj <- readRDS(file=paste0(opt_data$results_path,"s_obj_w_metcells.rds"))
    toc()
    # get list of each cell type
    celltype_list <- as.character(unique(seurat_obj@meta.data[[opt_data$level]]))
    n_metacells_per_group <- read.csv(paste0(opt_data$results_path, "n_metacells_per_group.csv"),row.names = 1)
    if (opt_filter_genes=="custom"){
      tic("loading network_gene_list.RData")
      network_gene_list <- readRDS(paste0(opt_data$results_path, "network_gene_list.RData"))
      toc()
    }
    toc()
    
    #loading s_obj_w_metcells.rds: 112.83 sec elapsed
    #loading network_gene_list.RData: 0.13 sec elapsed
    #loading metacells total: 112.96 sec elapsed
  } else{ #create and save seurat_obj with metacells
    tic("calculating and saving metacells total")
    seurat_obj <- get_seurat_object_ready(settings,opt_data)
    
    #some preparation before metacells can be generated:
    #to do: does this function work on the right assay and slot if sct norm data are used? 
    seurat_obj <- metacell_preparations(seurat_obj,settings,opt_data)
    
    # plot UMAP with annotated cell types
    # this is the same plot as module_features_hMEs and module_features_scores
    visualize_cell_type_UMAP(seurat_obj, opt_data$level, opt_data$results_path)
    
    #determine cell type markers if not already done  (checks if cell_type_markers.csv exists):
    markers <- define_cell_type_markers(seurat_obj,opt_data)
    
    # get list of each cell type
    celltype_list <- as.character(unique(seurat_obj@meta.data[[opt_data$level]]))
    
    if (opt_filter_genes=="custom"){
      # get list of genes for each cell type
      network_gene_list <- select_genes_for_network_analysis(seurat_obj, fraction = settings$fraction_of_cells_expr_gene, group.by = opt_data$level)
      tic("saving network_gene_list.RData")
      saveRDS(network_gene_list, file=paste0(opt_data$results_path, "network_gene_list.RData"))
      toc()
      
      all_network_genes <- unique(unlist(network_gene_list))
      
      display_genes_filtering_stats(celltype_list,all_network_genes,network_gene_list,seurat_obj,settings$fraction_of_cells_expr_gene)
      
      # only run metacell aggregation one time
      print('Setting up and running metacells.')
      
      #determine which slots in seurat_obj are used:
      seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = opt_filter_genes,
        gene_list = all_network_genes,#network_gene_list[[ct]],#shouldn't we have genes here expressed in at least one cell type and after metacell construction run SetupForWGCNA() for this Cell types genes?
        wgcna_name = "all_CTs_all_network_genes")
    } else{
      # only run metacell aggregation one time
      print('Setting up and running metacells.')
      
      #determine which slots in seurat_obj are used:
      seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = opt_filter_genes, # upcoming: explore other settings further
        wgcna_name = "all_CTs_all_network_genes")
    }
    
    #generate metacells by grouping by Sample and cell_type to achieve the desired result.
    #exclude underrepresented cell types with min_cells
    seurat_obj <- MetacellsByGroups(
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
      target_metacells = settings$target_number_meta_cells)
    
    #normalize metacell expression matrix:
    if (opt_data$normalize == "basic"){
      seurat_obj <- NormalizeMetacells(seurat_obj)
    } else if (opt_data$normalize == "metacells"){
      metacell_obj <- GetMetacellObject(seurat_obj)
      meta_data <- seurat_obj@meta.data[,-1] # gets everything that is stored under metadata except orig.ident
      metacell_obj <- normalize_data_with_sctransform(metacell_obj, opt_data)
      # only keep genes that were used for SCTransform
      
      if (packageDescription("Seurat")$Version >= package_version("5.0")){
        sct_data <- LayerData(metacell_obj, slot='scale.data', assay='SCT')
      } else{
        sct_data <- GetAssayData(metacell_obj, slot='scale.data', assay='SCT')
      }
      sct_genes <- rownames(sct_data)
      gene_list <- GetWGCNAGenes(seurat_obj)
      gene_list <- gene_list[gene_list %in% sct_genes]
      # update the genes used for WGCNA, and reset the metacell object
      seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list)
      seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj)
    }
    
    #visualize metacells:
    visualize_meta_cells(seurat_obj, settings, opt_data)
    # save information on metacells
    metacell_info <- seurat_obj@misc[["all_CTs_all_network_genes"]][["wgcna_metacell_obj"]]@meta.data
    #metacell_info <- seurat_obj@misc[[celltype_list[1]]][["wgcna_metacell_obj"]]@meta.data
    n_metacells_per_group <- as.data.frame.matrix(table(metacell_info[[opt_data$level]], metacell_info[[opt_data$group]]))
    write.csv(n_metacells_per_group, paste0(opt_data$results_path, "n_metacells_per_group.csv"))
    n_metacells_per_donor <- as.data.frame.matrix(table(metacell_info[[opt_data$level]], metacell_info[[opt_data$sample_variable]]))
    write.csv(n_metacells_per_donor, paste0(opt_data$results_path, "n_metacells_per_donor.csv"))
    #To Do: optionally: save seurat_obj with metacells
    if (opt_save_seurat_obj_with_metacells==TRUE){
      tic("saving s_obj_w_metcells.rds")
      saveRDS(seurat_obj, file=paste0(opt_data$results_path,"s_obj_w_metcells.rds"))
      toc()
      #saving s_obj_w_metcells.rds: 1053 sec elapsed
    }
    toc()
  } 
  return(list(seurat_obj,celltype_list,n_metacells_per_group,network_gene_list))
}

### Function to perform network analysis is now 3 separate functions
network_analysis <- function(seurat_obj,ct,ct_,opt_data,settings){
  
  #set up expression matrix:
  #for normalization == "counts" make sure the right assay and slot is used
  if (opt_data$normalize == "counts"){
    #sctransformed data is stored in a different assay:
    use_assay <- "SCT"
  } else{
    use_assay <- "RNA"
  }
  #slot is always data
  seurat_obj <- SetDatExpr_v2(
    seurat_obj,
    group_name = ct, # the name of the group of interest in the group.by column
    group.by=opt_data$level, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = use_assay, 
    slot = 'data', # using normalized data (NormalizeData() and scaleData() result instead of raw counts or sctransformed data instead of raw counts)
    use_metacells = TRUE #FALSE # it still uses the metacell info for building the network but extracts counts and metadata from RNA array, was the only way to make it run
  )
  
  #select soft-power threshold:
  # Test different soft powers:
  #hdWGCNA constructs a gene-gene correlation adjacency matrix to infer co-expression relationships between genes. 
  #The correlations are raised to a power to reduce the amount of noise present in the correlation matrix, thereby 
  #retaining the strong connections and removing the weak connections. Therefore, it is critical to determine a 
  #proper value for the soft power threshold.
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )
  #browser()
  visualize_soft_powers_results(seurat_obj,opt_data$results_path_ct_norm)
  
  power_table <- GetPowerTable(seurat_obj)
  write.csv(power_table, paste0(opt_data$results_path_ct_norm, "power_table.csv"))
  cat("power_table:\n")
  head(power_table)
  
  #construct co-expression network for current cell type:
  seurat_obj <- ConstructNetwork(
    seurat_obj, 
    #soft_power=1, #automatically detected by default
    tom_name = paste0(ct_,"__norm_",opt_data$normalize), # name of the topoligical overlap matrix written to disk
    wgcna_name = ct,
    setDatExpr=FALSE,
    randomSeed = settings$randomSeed,
    overwrite_tom = TRUE,
    TomType = settings$TOMType,#
    networkType = settings$networkType,
    mergeCutHeight = settings$mergeCutHeight,
    corType = settings$corType,
    detectCutHeight = settings$detectCutHeight,
    deepSplit = settings$deepSplit,
    minModuleSize = settings$minModuleSize
  )
  
  
  #plot the dendrogram
  visualize_dendrogram(seurat_obj,opt_data$results_path_ct_norm,ct)
  
  #inspect the topological overlap matrix
  #path was somehow setup in the wrong way, correct this mistake:
  if (settings$export_TOM_matrix==TRUE){
    #TOM_path_with_error <- seurat_obj@misc$scz$wgcna_net$TOMFiles
    #filename = str_split(TOM_path_with_error,'//',n=2)
    #seurat_obj@misc$scz$wgcna_net$TOMFiles <- paste0(opt_data$results_path, filename[[1]][2])
    TOM <- GetTOM(seurat_obj, wgcna_name = ct) # gene by gene matrix
    if (opt_data$opt_exclude_MT_genes==TRUE){
      res_path <- paste0(opt_data$code_path,"TOM/without_MT_genes/")
    } else{
      res_path <- paste0(opt_data$code_path,"TOM/")
    }
    #make sure res_path exists:
    if (!dir.exists(res_path)){
      dir.create(res_path)
    } 
    saveRDS(TOM, file=paste0(res_path,ct_,"__norm_",opt_data$normalize,"_TOM.Rds"))
    
    #write.table(TOM, file=paste0(ct_,"__norm_",opt_data$normalize,"_TOM.Rdata"))
    modules <- GetModules(seurat_obj) %>% 
                          subset(module != 'grey') %>% 
                          mutate(module = droplevels(module))
    #write.table(modules, file=paste0(ct_,"__norm_",opt_data$normalize,"_modules_for_TOM.Rdata"))
    if (opt_data$opt_exclude_MT_genes==TRUE){
      saveRDS(modules, file=paste0(opt_data$code_path,"TOM/without_MT_genes/",ct_,"__norm_",opt_data$normalize,"_modules_for_TOM.Rds"))
    } else{
      saveRDS(modules, file=paste0(opt_data$code_path,"TOM/",ct_,"__norm_",opt_data$normalize,"_modules_for_TOM.Rds"))
    }
  }
  
  return(seurat_obj)
}

network_analysis_exploration <- function(seurat_obj,ct,ct_,opt_data,settings,exploration_settings_name){
  #browser()
  #set up expression matrix:
  #for normalization == "counts" make sure the right assay and slot is used
  if (opt_data$normalize == "counts"){
    #sctransformed data is stored in a different assay:
    use_assay <- "SCT"
  } else{
    use_assay <- "RNA"
  }
  #slot is always data
  seurat_obj <- tryCatch(
    {
      seurat_obj <- SetDatExpr_v2(
        seurat_obj,
        group_name = ct, # the name of the group of interest in the group.by column
        group.by=opt_data$level, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
        assay = use_assay, 
        slot = 'data', # using normalized data (NormalizeData() and scaleData() result instead of raw counts or sctransformed data instead of raw counts)
        use_metacells = TRUE #FALSE # it still uses the metacell info for building the network but extracts counts and metadata from RNA array, was the only way to make it run
      )
    },
    error=function(cond){
      message("SetDatExpr failed!")
      message(cond)
      return(c())
    }
  )
  if (length(seurat_obj)==0){
    return(c())
  }
  
  #select soft-power threshold:
  # Test different soft powers:
  #hdWGCNA constructs a gene-gene correlation adjacency matrix to infer co-expression relationships between genes. 
  #The correlations are raised to a power to reduce the amount of noise present in the correlation matrix, thereby 
  #retaining the strong connections and removing the weak connections. Therefore, it is critical to determine a 
  #proper value for the soft power threshold.
  #browser()
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )
  
  #browser()
  #construct co-expression network for current cell type:
  
  seurat_obj <- tryCatch(
    {
      seurat_obj <- ConstructNetwork(
        seurat_obj, 
        tom_name = paste0(ct_,"__norm_",opt_data$normalize,'_',exploration_settings_name), # name of the topoligical overlap matrix written to disk
        wgcna_name = exploration_settings_name,
        metacell_location = "all_CTs_all_network_genes",
        #setDatExpr=FALSE, # what does this do? can't find it mentioned in any function
        randomSeed = settings$randomSeed,
        overwrite_tom = TRUE,
        TomType = settings$TOMType,#
        networkType = settings$networkType,
        mergeCutHeight = settings$mergeCutHeight,
        corType = settings$corType,
        detectCutHeight = settings$detectCutHeight,
        deepSplit = settings$deepSplit,
        minModuleSize = settings$minModuleSize
      )
    },
    error=function(cond){
      message("ConstructNetwork failed!")
      message(cond)
      return(c())
    }
  )
  if (length(seurat_obj)==0){
    return(c())
  }
  
  #plot the dendrogram
  graphic_filename <- paste0(opt_data$results_path_long,"dendrogram.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  PlotDendrogram(seurat_obj, main=exploration_settings_name)
  dev.off()

  return(seurat_obj)
}

module_analyis <- function(seurat_obj,opt_data,ct,settings){
  #Module Eigengenes and Connectivity
  #commonly used metric to summarize the gene expression profile of an entire co-expression module
  # module eigengenes are computed by performing principal component analysis (PCA) on the subset 
  #of the gene expression matrix comprising each module. The first PC of each of these PCA matrices 
  #are the MEs.
  
  # need to run ScaleData first or else harmony throws an error:
  #seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  
  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj, harmonized=TRUE)
  write.csv(hMEs, paste0(opt_data$results_path_ct_norm, "hMEs.csv"))
  
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  write.csv(MEs, paste0(opt_data$results_path_ct_norm, "MEs.csv"))
  
  #module connectivity:
  #want to focus on the "hub genes", those which are highly connected within each module.
  #eigengene-based connectivity: kME of each gene
  #ModuleConnectivity to compute the kME values in the full single-cell dataset, rather than the metacell dataset
  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = opt_data$level, group_name = ct
  )
  
  # rename the modules
  # seurat_obj <- ResetModuleNames(
  #   seurat_obj,
  #   new_name = "M"
  # )
  
  # plot genes ranked by kME for each module
  visualize_modules(seurat_obj,opt_data$results_path_ct_norm)
  
  #Getting the module assignment table:
  # get the module assignment table:
  modules <- GetModules(seurat_obj)
  write.csv(modules, paste0(opt_data$results_path_ct_norm, "modules_df.csv") )
  
  # show the first up to 6 columns or all:
  cat("modules:\n")
  head(modules[,1:min(dim(modules)[2],6)])
  
  # table of the top N hub genes sorted by kME:
  hub_df <- GetHubGenes_v2(seurat_obj, n_hubs = 10)
  write.csv(hub_df, paste0(opt_data$results_path_ct_norm, "hub_df.csv"))
  hub_df <- GetHubGenes_v2(seurat_obj, n_hubs = "all")
  write.csv(hub_df, paste0(opt_data$results_path_ct_norm, "hub_df_all.csv"))
  
  #Compute hub gene signature scores
  #Gene scoring analysis: computing a score for the overall signature of a set of genes
  #function ModuleExprScore computes gene scores for a give number of genes for each module, 
  #using either the Seurat or UCell algorithm. 
  #Gene scoring is an alternative way of summarizing the expression of a module from computing the module eigengene.
  # compute gene scoring for the top 25 hub genes by kME for each module
  # with Seurat method
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method=settings$gene_sign_score_comp 
  )
  
  #Module feature plot: 
  visualize_module_features(seurat_obj,opt_data$results_path_ct_norm, features_to_visualize = 'hMEs')
  visualize_module_features(seurat_obj,opt_data$results_path_ct_norm, features_to_visualize = 'scores')
  
  # plot module correlagram
  seurat_obj <- AvgModuleExpr(
    seurat_obj, 
    n_genes = 25, 
    wgcna_name = ct
  )
  visualize_module_correlations(seurat_obj,opt_data$results_path_ct_norm,ct)
  
  #add hMEs to seurat meta_data
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)
  mods <- colnames(hMEs)
  mods <- mods[mods != 'grey']
  
  visualize_hMEs_as_dotplot(seurat_obj,mods,opt_data)
  
  # To DO: use modulenames stored in seurat_obj for loop
  for (module_name in mods) {
    visualize_hMEs_as_violins(seurat_obj,opt_data$level,opt_data$results_path_ct_norm, module_name)
  }
  visualize_module_gene_expression_as_violins(seurat_obj,opt_data$level,opt_data$results_path_ct_norm, modules)
  
  # module trait correlation analysis
  # convert character traits to factors
  for (trait in settings$traits) {
    if (is.character(seurat_obj@meta.data[[trait]])) {
      seurat_obj[[trait]] <- lapply(seurat_obj[[trait]], factor) 
    }
  }
  
  seurat_obj <- ModuleTraitCorrelation_v2(
    seurat_obj,
    traits = settings$traits,
    group.by=opt_data$level,
    subset_by = opt_data$level,
    subset_groups = celltype_list # only include other cell types included in network analysis
  )
  
  visualize_module_trait_corr(seurat_obj,opt_data$results_path_ct_norm)
  
  # Marker gene overlap analysis
  markers <- read.csv(paste0(opt_data$results_path, "cell_type_markers.csv"))
  overlap_df <- OverlapModulesDEGs(
    seurat_obj,
    deg_df = markers,
    fc_cutoff = 1 # log fold change cutoff for overlap analysis
  )
  # save overlap with cell type markers
  write.csv(overlap_df, row.names=FALSE, quote=FALSE, file=paste0(opt_data$results_path_ct_norm, "module_overlap_with_cell_type_markers.csv"))
  # visualize overlap
  visualize_marker_overlap(overlap_df,opt_data$results_path_ct_norm,ct)
  
  return(seurat_obj)
}

run_gprofiler_for_network_modules <- function(seurat_obj,ct_,opt_data){
  
  # load background for gprofiler (here: all genes expressed in our dataset)
  g_profiler_background <- read.csv(opt_data$data_background_path, header = TRUE)$ensgid
  
  # get translator between gene names and accessions for later 
  translator <- get_translator(seurat_obj)
  
  #get modules:
  modules <- GetModules(seurat_obj)
  hMEs <- GetMEs(seurat_obj)
  
  module_gene_mapping <- save_module_gene_mapping(module_detect=modules, cell_type=ct_, translator=translator, results_dir = opt_data$results_path)
  
  for (module_name in colnames(hMEs)) {
    #run gprofiler and store results:
    run_gprofiler(module_gene_mapping, module_name, opt_data$results_path_ct_norm, ct_, g_profiler_background)
  }
}

save_module_gene_mapping <- function(module_detect. = module_detect, cell_type. = cell_type, translator. = translator, results_dir. = results_dir){
  # prepare for saving some more infos
  #if used by hdWGCNA_main: module_detect.$color
  #previously for WGCNA_main: module_detect.$colors
  module_gene_mapping <- as.data.frame(module_detect.$color, row.names=module_detect.$gene_name) %>%
    mutate(module_name = paste0(cell_type., "_module_", as.character(`module_detect.$color`))) %>%
    select(-`module_detect.$color`)
  module_gene_mapping$gene_ids <- rownames(module_gene_mapping)
  module_gene_mapping <- inner_join(module_gene_mapping, translator., by="gene_ids")
  
  write_csv(module_gene_mapping, paste0(results_dir., "module_gene_mapping_info_", cell_type., ".csv"))
  
  return(module_gene_mapping)
}

run_gprofiler <- function(module_gene_mapping, module, results_dir, cell_type,g_profiler_background){
  module_list <- c()
  for (module in unique(module_gene_mapping$module_name)) {
    if (module!=paste0(cell_type, "_module_", as.character(0))) { # grey not included since it includes all the "rest" genes that are not part of a different module
      module_list <- c(module_list, paste0("> ", module))
      module_list <- c(module_list, module_gene_mapping[module_gene_mapping$module_name==module,]$accession)
    }
  }
  g_profiler_df <- as.data.frame(module_list)
  write.table(g_profiler_df, file = paste0(results_dir, "module_file_for_gprofiler_", cell_type, ".csv"), quote=FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
  
  # perform enrichment analysis with g profiler here
  g_profiler_list <- list()
  for (module in sort(unique(module_gene_mapping$module_name))) {
    if (module!=paste0(cell_type, "_module_", as.character(0))) { # grey not included since it includes all the "rest" genes that are not part of a different module
      g_profiler_list[[module]] <- c(module_gene_mapping[module_gene_mapping$module_name == module,]$accession)
    }
  }

  organism <- "hsapiens"

  # TF removed from source list
  source_list <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP")
  gostres <- gost(query = g_profiler_list, 
                  organism =organism, ordered_query = FALSE, # is it a ranked gene list?
                  multi_query = TRUE, significant = TRUE, # returns only significant results
                  exclude_iea = TRUE, # exclude electronic GO annotations?#was previously FALSE
                  measure_underrepresentation = FALSE, #measure under-representation in gene set?
                  evcodes = FALSE, 
                  user_threshold = 0.05, 
                  correction_method = "g_SCS", #was previously "bonferroni", 
                  domain_scope = "custom_annotated", # was previously "annotated", 
                  custom_bg = g_profiler_background, # custom background data set - change this to only include genes we detect (more brain focused)
                  numeric_ns = "", 
                  as_short_link = FALSE,
                  sources = source_list)
  
  GO_results <- gostres$result
  GO_results <- GO_results[ , -which(names(GO_results) == "parents")] # remove this to be able to save as csv file
  
  # split up the p_value column
  # WHY "P_VALUE"? IS IT A BONFERRONI CORRECTED Q_VALUE OR NOT? --> Have to double-check!
  p_val_df <- t(as.data.frame(GO_results$p_values))
  colnames(p_val_df) <- paste0("pval_", names(g_profiler_list))                      # Order should be the same since g_profiler_list was sorted when created 
  rownames(p_val_df) <- c()
  GO_results <- cbind(GO_results, p_val_df)
  
  # save gprofiler results as csv 
  write_csv(GO_results, paste0(results_dir, 'g_profiler_results_', cell_type, '.csv'))
  
  # get top 10 hits by pvalue for each module
  # right now, the category is ignored, maybe the top 3 for each module-category-combination is interesting too
  highlights <- c()
  for (i in paste0("pval_", names(g_profiler_list))) {
    eval(parse(text = paste0("GO_results_sel <- GO_results[order(GO_results$\'", i, "\'),]")))
    highlights <- c(highlights, GO_results_sel$term_id[1:10])
  }
  highlights <- unique(highlights)
  
  # get gprofiler results image and save as pdf 
  p = gostplot(gostres, capped = FALSE, interactive = FALSE)
  publish_gostplot(p = p,
                   highlight_terms = highlights,
                   filename = paste0(results_dir, 'g_profiler_results_', cell_type, '.pdf'),
                   width = 20, height = 30)
}

create_ct_subfolder <- function(opt_data,ct_){
  results_path_ct_norm <- paste0(opt_data$results_path,'/',ct_,'/',"norm_",opt_data$normalize,"/")
  if (file.exists(results_path_ct_norm)==FALSE) {
    dir.create(results_path_ct_norm, recursive=TRUE)
  }
  return(results_path_ct_norm)
}

### Functions used by metacells_settings_exploration_main.R
# Function to create folder for results 
create_subfolder <- function(opt_data, k,s,m){
  subfolder_str <- paste0(as.character(k),'_neighbors_',as.character(s),'_cells_shared_at_least_',as.character(m),'_cells/')
  opt_data$results_path_long <- paste0(opt_data$results_path,subfolder_str)
  if (file.exists(opt_data$results_path_long)==FALSE) {
    dir.create(opt_data$results_path_long, recursive=TRUE)
  }
  return(opt_data)
}

#Function used by dendrogram settings_exploration_main.R
create_module_subfolder <- function(opt_data,k,s,m){
  #subfolder_str <- paste0(as.character(k),'_deepSplit_',as.character(s),'_mergeCutHeight_',as.character(m),'_minModuleSize/')
  subfolder_str <- paste0('deepSplit_',as.character(k),'_mergeCutHeight_', gsub('.', '', as.character(s), fixed = TRUE),'_minModuleSize_',as.character(m),'/')
  opt_data$results_path_long <- paste0(opt_data$results_path_ct_norm,subfolder_str)
  if (file.exists(opt_data$results_path_long)==FALSE) {
    dir.create(opt_data$results_path_long, recursive=TRUE)
  }
  return(opt_data)
}

# Functions to assess success of metacell calculation
calculate_sparsity <- function(M){
  SPARSITY <- sum(M==0)/(dim(M)[1]*dim(M)[2])
  return(SPARSITY)
}

#taken from from https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r:
sparse.cor3 <- function(x){
  #memory.limit(size=10000)
  n <- nrow(x)
  
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  C <- covmat/crossprod(t(sdvec)) # correlation matrix
  
  return(C)
}

### Functions adopted from Paul's WGCNA code

# get translator between gene names and accessions for later 
get_translator <- function(data.seurat){
  translator_raw <- data.seurat@assays[["RNA"]]@meta.features
  translator <- data.frame(cbind(translator_raw$Accession, rownames(translator_raw)))
  colnames(translator) <- c('accession','gene_ids')
  
  return(translator)
}

DME_analysis <- function(seurat_obj, ct, ct_, opt_data, settings, DME_list, opt_test) {
  case_group <- 'SCZ'

  # case group: SCZ or mutant
  group1 <- seurat_obj@meta.data %>% subset(seurat_obj@meta.data[opt_data$level] == ct & seurat_obj@meta.data[opt_data$group] == case_group) %>% rownames
  # control group
  group2 <- seurat_obj@meta.data %>% subset(seurat_obj@meta.data[opt_data$level] == ct & seurat_obj@meta.data[opt_data$group] != case_group) %>% rownames
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1, # positive fold-change means up-regulated in this group
    barcodes2 = group2, # negative fold-change means up-regulated in this group
    test.use=opt_test,
    wgcna_name=ct
  )
  
  
  # plot figures
  visualize_DMEs_as_lollipop(seurat_obj, DMEs, settings, ct, ct_, opt_data$results_path_ct_norm, opt_test)
  visualize_DMEs_as_volcano(seurat_obj, DMEs, settings, ct, ct_, opt_data$results_path_ct_norm, opt_test)
  
  DMEs$level <- ct
  
  # fix infs:
  DMEs$avg_log2FC <- ifelse(abs(DMEs$avg_log2FC) == Inf, 0, DMEs$avg_log2FC)
  DME_list[[ct]] <- DMEs
  
  return(DME_list)
}

ModuleTraitCorrelation_v2 <- function (seurat_obj, traits, group.by = NULL, features = "hMEs", 
                                       cor_method = "pearson", subset_by = NULL, subset_groups = NULL, 
                                       wgcna_name = NULL, ...) 
{
  #browser()
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (features == "hMEs") {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  }
  else if (features == "MEs") {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  }
  else if (features == "scores") {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  }
  else {
    stop("Invalid feature selection. Valid choices: hMEs, MEs, scores, average")
  }
  if (!is.null(subset_by)) {
    print("subsetting")
    seurat_full <- seurat_obj
    MEs <- MEs[seurat_obj@meta.data[[subset_by]] %in% subset_groups, 
    ]
    seurat_obj <- seurat_obj[, seurat_obj@meta.data[[subset_by]] %in% 
                               subset_groups]
  }
  if (sum(traits %in% colnames(seurat_obj@meta.data)) != length(traits)) {
    stop(paste("Some of the provided traits were not found in the Seurat obj:", 
               paste(traits[!(traits %in% colnames(seurat_obj@meta.data))], 
                     collapse = ", ")))
  }
  if (is.null(group.by)) {
    group.by <- "temp_ident"
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }
  valid_types <- c("numeric", "factor", "integer")
  data_types <- sapply(traits, function(x) {
    class(seurat_obj@meta.data[, x])
  })
  if (!all(data_types %in% valid_types)) {
    incorrect <- traits[!(data_types %in% valid_types)]
    stop(paste0("Invalid data types for ", paste(incorrect, 
                                                 collapse = ", "), ". Accepted data types are numeric, factor, integer."))
  }
  if (any(data_types == "factor")) {
    factor_traits <- traits[data_types == "factor"]
    for (tr in factor_traits) {
      warning(paste0("Trait ", tr, " is a factor with levels ", 
                     paste0(levels(seurat_obj@meta.data[, tr]), collapse = ", "), 
                     ". Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?"))
    }
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  trait_df <- seurat_obj@meta.data[, traits]
  if (length(traits == 1)) {
    trait_df <- data.frame(x = trait_df)
    colnames(trait_df) <- traits
  }
  if (any(data_types == "factor")) {
    factor_traits <- traits[data_types == "factor"]
    for (tr in factor_traits) {
      trait_df[, tr] <- as.numeric(trait_df[, tr])
    }
  }
  cor_list <- list()
  pval_list <- list()
  fdr_list <- list()
  temp <- Hmisc::rcorr(as.matrix(trait_df), as.matrix(MEs), 
                       type = cor_method)
  cur_cor <- temp$r[traits, mods]
  cur_p <- temp$P[traits, mods]
  p_df <- cur_p %>% reshape2::melt()
  if (length(traits) == 1) {
    tmp <- rep(mods, length(traits))
    tmp <- factor(tmp, levels = mods)
    tmp <- tmp[order(tmp)]
    p_df$Var1 <- traits
    p_df$Var2 <- tmp
    rownames(p_df) <- 1:nrow(p_df)
    p_df <- dplyr::select(p_df, c(Var1, Var2, value))
  }
  p_df <- p_df %>% dplyr::mutate(fdr = p.adjust(value, method = "fdr")) %>% 
    dplyr::select(c(Var1, Var2, fdr))
  cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = "fdr")
  rownames(cur_fdr) <- cur_fdr$Var1
  cur_fdr <- cur_fdr[, -1]
  cor_list[["all_cells"]] <- cur_cor
  pval_list[["all_cells"]] <- cur_p
  fdr_list[["all_cells"]] <- cur_fdr
  trait_df <- cbind(trait_df, seurat_obj@meta.data[, group.by])
  colnames(trait_df)[ncol(trait_df)] <- "group"
  MEs <- cbind(as.data.frame(MEs), seurat_obj@meta.data[, 
                                                        group.by])
  colnames(MEs)[ncol(MEs)] <- "group"
  if (class(seurat_obj@meta.data[, group.by]) == "factor") {
    group_names <- levels(seurat_obj@meta.data[, group.by])
  }
  else {
    group_names <- levels(as.factor(seurat_obj@meta.data[, 
                                                         group.by]))
  }
  trait_list <- dplyr::group_split(trait_df, group, .keep = FALSE)
  ME_list <- dplyr::group_split(MEs, group, .keep = FALSE)
  names(trait_list) <- group_names
  names(ME_list) <- group_names
  for (i in names(trait_list)) {
    temp <- Hmisc::rcorr(as.matrix(trait_list[[i]]), as.matrix(ME_list[[i]]))
    cur_cor <- temp$r[traits, mods]
    cur_p <- temp$P[traits, mods]
    p_df <- cur_p %>% reshape2::melt()
    if (length(traits) == 1) {
      tmp <- rep(mods, length(traits))
      tmp <- factor(tmp, levels = mods)
      tmp <- tmp[order(tmp)]
      p_df$Var1 <- traits
      p_df$Var2 <- tmp
      rownames(p_df) <- 1:nrow(p_df)
      p_df <- dplyr::select(p_df, c(Var1, Var2, value))
    }
    p_df <- p_df %>% dplyr::mutate(fdr = p.adjust(value, 
                                                  method = "fdr")) %>% dplyr::select(c(Var1, Var2, 
                                                                                       fdr))
    cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = "fdr")
    rownames(cur_fdr) <- cur_fdr$Var1
    cur_fdr <- cur_fdr[, -1]
    cor_list[[i]] <- cur_cor
    pval_list[[i]] <- cur_p
    fdr_list[[i]] <- as.matrix(cur_fdr)
  }
  mt_cor <- list(cor = cor_list, pval = pval_list, fdr = fdr_list)
  if (!is.null(subset_by)) {
    seurat_full <- SetModuleTraitCorrelation(seurat_full, 
                                             mt_cor, wgcna_name)
    seurat_obj <- seurat_full
  }
  else {
    seurat_obj <- SetModuleTraitCorrelation(seurat_obj, 
                                            mt_cor, wgcna_name)
  }
  seurat_obj
}

display_genes_filtering_stats <- function(celltype_list,all_network_genes,network_gene_list,seurat_obj,fraction){
  print(paste0("number of genes for fraction = ",as.character(fraction),": "))
  for (ct in celltype_list){
    print(paste0("for ",ct,": ",as.character(eval(parse(text=paste0("length(network_gene_list$`",ct,"`)"))))))
  }
  n_genes_kept <- length(all_network_genes)
  print(paste0("for all cell types: ",n_genes_kept))
  n_genes_rem <- length(seurat_obj@assays$RNA@counts@Dimnames$features) - n_genes_kept
  p_genes_rem <- round(100*(n_genes_rem/length(seurat_obj@assays$RNA@counts@Dimnames$features)))
  print(paste0("genes not used for metacells: ",n_genes_rem," (~",p_genes_rem,"%)"))
}

#save seurat object:
save_seurat_obj <- function(seurat_obj,opt_data,ct_){
  #make sure path exists:
  if (opt_data$opt_exclude_MT_genes == TRUE){
    path_seuratObj = paste0(opt_data$code_path,"/seurat_objects/without_MT_genes/")
  } else{
    path_seuratObj = paste0(opt_data$code_path,"/seurat_objects/")
  }
  if (file.exists(path_seuratObj)==FALSE) {
    dir.create(path_seuratObj, recursive=TRUE)
  }
  SaveH5Seurat(object=seurat_obj,
               filename=paste0(path_seuratObj,"S_",ct_),
               overwrite=TRUE,
               verbose=TRUE)
}

#allow n_hubs to be "all" and export all genes in the network
GetHubGenes_v2 <- function(seurat_obj, n_hubs = 10, mods = NULL, wgcna_name = NULL) 
  {
    if (is.null(wgcna_name)) {
      wgcna_name <- seurat_obj@misc$active_wgcna
    }
    modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 
                                                               "grey")
    if (is.null(mods)) {
      mods <- levels(modules$module)
      mods <- mods[mods != "grey"]
    }
    else {
      if (!all(mods %in% modules$module)) {
        stop("Invalid selection for mods.")
      }
    }
    hub_df <- do.call(rbind, lapply(mods, function(cur_mod) {
                                            cur <- subset(modules, module == cur_mod)
                                            cur <- cur[, c("gene_name", "module", paste0("kME_", 
                                                                                         cur_mod))]
                                            names(cur)[3] <- "kME"
                                            cur <- dplyr::arrange(cur, kME)
                                            if (n_hubs=="all"){
                                              cur %>% dplyr::top_frac(n=1, wt=kME)
                                            } else{
                                              cur %>% dplyr::top_n(n_hubs, wt = kME)
                                            }
                                          }
                                    )
                      )
    if (nrow(hub_df) >=1){
      rownames(hub_df) <- 1:nrow(hub_df)
    }
    hub_df
  }

#getAssayData() replaced in seurat with LayerData()
SetDatExpr_v2 <- function (seurat_obj, group_name, use_metacells = TRUE, group.by = NULL, 
            multi.group.by = NULL, multi_group_name = NULL, return_seurat = TRUE, 
            assay = NULL, slot = "data", mat = NULL, wgcna_name = NULL, 
            ...) 
  {
    if (is.null(wgcna_name)) {
      wgcna_name <- seurat_obj@misc$active_wgcna
    }
    if (!CheckWGCNAName(seurat_obj, wgcna_name)) {
      stop(paste0("Invalid wgcna_name supplied: ", wgcna_name))
    }
    if (is.null(assay)) {
      assay <- DefaultAssay(seurat_obj)
      warning(paste0("assay not specified, trying to use assay ", 
                     assay))
    }
    if (!(assay %in% Assays(seurat_obj))) {
      stop(paste0("Invalid choice of assay: ", assay, " not found in Assays(seurat_obj)."))
    }
    if (!(slot %in% c("counts", "data", "scale.data"))) {
      stop("Invalid choice of slot. Valid choices are counts, data, or scale.data.")
    }
    params <- GetWGCNAParams(seurat_obj, wgcna_name)
    genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
    if (is.null(mat)) {
      m_obj <- GetMetacellObject(seurat_obj, wgcna_name)
      if (use_metacells & !is.null(m_obj)) {
        s_obj <- m_obj
      }
      else {
        if (is.null(m_obj)) {
          warning("Metacell Seurat object not found. Using full Seurat object instead.")
        }
        s_obj <- seurat_obj
      }
      seurat_meta <- s_obj@meta.data
      if (!is.null(group.by)) {
        if (!(group.by %in% colnames(s_obj@meta.data))) {
          m_cell_message <- ""
          if (use_metacells) {
            m_cell_message <- "metacell"
          }
          stop(paste0(group.by, " not found in the meta data of the ", 
                      m_cell_message, " Seurat object"))
        }
        if (!all(group_name %in% s_obj@meta.data[[group.by]])) {
          groups_not_found <- group_name[!(group_name %in% 
                                             s_obj@meta.data[[group.by]])]
          stop(paste0("Some groups in group_name are not found in the seurat_obj: ", 
                      paste(groups_not_found, collapse = ", ")))
        }
      }
      if (!is.null(multi.group.by)) {
        if (!(multi.group.by %in% colnames(s_obj@meta.data))) {
          m_cell_message <- ""
          if (use_metacells) {
            m_cell_message <- "metacell"
          }
          stop(paste0(multi.group.by, " not found in the meta data of the ", 
                      m_cell_message, " Seurat object"))
        }
        if (!all(multi_group_name %in% s_obj@meta.data[[multi.group.by]])) {
          groups_not_found <- multi_group_name[!(multi_group_name %in% 
                                                   s_obj@meta.data[[multi.group.by]])]
          stop(paste0("Some groups in group_name are not found in the seurat_obj: ", 
                      paste(groups_not_found, collapse = ", ")))
        }
      }
      if (!is.null(group.by)) {
        seurat_meta <- seurat_meta %>% subset(get(group.by) %in% 
                                                group_name)
      }
      if (!is.null(multi.group.by)) {
        seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in% 
                                                multi_group_name)
      }
      cells <- rownames(seurat_meta)
      
      if (packageDescription("Seurat")$Version >= package_version("5.0")){
        datExpr <- as.data.frame(Seurat::LayerData(s_obj, 
                                                   assay = assay, slot = slot)[genes_use, cells])
      } else{
        datExpr <- as.data.frame(Seurat::GetAssayData(s_obj, 
                                                      assay = assay, slot = slot)[genes_use, cells])
      }
      
      datExpr <- as.data.frame(t(datExpr))
    }
    else {
      datExpr <- mat
      if (any(class(datExpr) != "data.frame")) {
        datExpr <- as.data.frame(datExpr)
      }
      if (!all(colnames(datExpr) %in% rownames(seurat_obj))) {
        stop("colnames of the provided matrix are invalid. Make sure that the colnames are features (genes), and that all of these features are in the seurat_obj")
      }
      datExpr <- datExpr[, genes_use]
    }
    if (return_seurat) {
      gene_list <- genes_use[WGCNA::goodGenes(datExpr, ...)]
      datExpr <- datExpr[, gene_list]
      seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)
      seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
      out <- seurat_obj
    }
    else {
      out <- datExpr
    }
    out
  }

#GSA on modules
get_gene_lists_to_test <- function(path_genelist_to_test){
  CT_directories <-  list.dirs(path_genelist_to_test,recursive=FALSE)
  CT_directories[grep("/patient_clustering",CT_directories,value=FALSE)] <- ""
  CT_directories[grep("/figures",CT_directories,value=FALSE)] <- ""
  CT_directories[grep("/GSA_results",CT_directories,value=FALSE)] <- ""
  CT_directories <- CT_directories[CT_directories!=""]
  
  bool_first <- TRUE
  
  for (ct_folder in CT_directories){
    #determine which modules are significant
    M_sign <- read.csv(paste0(ct_folder,"/norm_basic/","DMEs_LR.csv"))
    ct <- unique(M_sign$level)
    ct <- str_replace_all(ct," ","_")
    Ms <- M_sign$module[M_sign$p_val_adj<0.05]
    M_genes <- read.csv(paste0(path_genelist_to_test,"module_gene_mapping_info_",ct,".csv"))
    for (m in Ms){
      mn <- unique(M_genes$module_name[endsWith(M_genes$module_name,m)])
      mn <- str_replace_all(mn,"-","_")
      #eval(parse(text=paste0(mn ,"<- M_genes$accession[endsWith(M_genes$module_name,m)]")))
      gene_list_to_test <- M_genes[endsWith(M_genes$module_name,m),]
      #add column of bools
      eval(parse(text=paste0("gene_list_to_test$",mn," <- TRUE")))
      #remove columns that are not needed
      gene_list_to_test <- subset(gene_list_to_test,select=-module_name)
      gene_list_to_test <- subset(gene_list_to_test,select=-gene_ids)
      if (bool_first){
        gene_list_to_test_bool <- gene_list_to_test
        bool_first <- FALSE
      } else{
        gene_list_to_test_bool <- merge(gene_list_to_test_bool,gene_list_to_test,by="accession",all.x = TRUE, all.y=TRUE)
      }
    }
  }
  
  #make accessions rownames
  row.names(gene_list_to_test_bool) = gene_list_to_test_bool$accession
  
  #rename accession to ensgid
  colnames(gene_list_to_test_bool)[colnames(gene_list_to_test_bool)=="accession"] <- "ensgid"
  
  #replace NA with FALSE:
  for (mn in colnames(gene_list_to_test_bool)){
    eval(parse(text=paste0("gene_list_to_test_bool$",mn,"[is.na(gene_list_to_test_bool$",mn,")] <- FALSE")))
  }
  
  return(gene_list_to_test_bool)
}
