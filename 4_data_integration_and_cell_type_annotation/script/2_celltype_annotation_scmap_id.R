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
main_path = getwd()
setwd("../")
main_project_path = getwd()

## settings:
opt_server = FALSE
opt_test = FALSE
alignment_method = 'cellranger'
opt_load_from_loom = TRUE #FALSE
opt_save_as_loom = TRUE# FALSE
opt_doublet_version = '_'
#should high drop out rate genes be specified?
opt_determine_HDR_genes=FALSE

if (opt_test==TRUE){
  opt_scmap <- c('cluster','cell2cluster','cell')
} else{
  opt_scmap <- 'cell2cluster'#'cluster' #
}

ref_data = "Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq"
ref_data_folder = "ABM_MCA/"
ref_data_subfolder_str="red_to_"
n_CT = 51 #18 or 51 or 76 or 120

cores=detectCores() # setup parallel backend to use many processors
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

if (opt_test==TRUE){
  n_cells_query <- 2500 # number of cells in test data set
  n_genes_ref <- c(500,1000,2000,3000,5000) # number of genes selected for reference data set
  red_method <- c('HDO','HVG_with_cell_ranger','HVG_with_seurat_v3')# method used for reducing number of genes in reference data set
} else {
  n_cells_query <- c(2500) # number of cells in query files
  n_genes_ref <- c(5000) #c(2000) #number of genes selected for reference data set
  red_method <- c('HVG_with_cell_ranger') #method used for reducing number of genes in reference data set
}

## 1) specify paths and load function to read data
path_code <- paste0(main_path_project,'/4_data_integration_and_cell_type_annotation/script/')

if (opt_test==TRUE){
  path_query_data <- paste0(main_path_project,'/4_data_integration_and_cell_type_annotation/data/reference_data_sets/',ref_data_folder)
  path_results <- paste0(main_path_project,'/4_data_integration_and_cell_type_annotation/output/CT_anno_scmap/Test_',ref_data_folder)
} else {
  path_query_data <- paste0(main_path_project,'/4_data_integration_and_cell_type_annotation/output/')
  path_results <- paste0(main_path_project,'/4_data_integration_and_cell_type_annotation/output/CT_anno_scmap/Filtered_',ref_data_folder,alignment_method,'/')
}

path_results_figures <- paste0(path_results,n_CT,"_CTs/",'figures',opt_doublet_version,'/')
path_results_files <- paste0(path_results,n_CT,"_CTs/",'output_files',opt_doublet_version,'/')

## 2) specify file names
if (opt_test==TRUE){
  query_filelist <- paste('matrix_',n_cells_query,sep="",'_cells_test.csv',sep="")
  query_file_GT <- paste('matrix_',n_cells_query,sep="",'_cells_GT_test.csv',sep="")
} else {
  pattern_query = paste0("S.*_",alignment_method,"_TH_and_D_adj",opt_doublet_version,"_filtered")
  if (opt_load_from_loom==TRUE){
    query_filelist <- list.files(path_query_data, pattern=paste0(pattern_query,".loom"), all.files=FALSE,full.names=FALSE)
  } else{
    query_filelist <- list.files(path_query_data, pattern=paste0(pattern_query,"_.*.csv"), all.files=FALSE,full.names=FALSE)
  }
}
ref_file_symbol_name <- '9606_symbol.csv'
ref_file_map_name <- paste(paste('9606_map_',n_CT,sep=""),'.csv',sep="")

## 3) initilize output variables in case of test option
if (opt_test==TRUE){
  settings <- character(0)
  kappa_val <- numeric(0)
  p_unassigned <- numeric(0)
}

subpath_ref_data <- paste0(main_path,'4_data_integration_and_cell_type_annotation/output/reference_data_sets/',ref_data_folder,ref_data_subfolder_str,as.character(n_CT))

## 4) annotate cell types
for (m in red_method){
  for (ng in n_genes_ref){
    path_ref_data <- paste0(subpath_ref_data,ng,'_',m,'/')
    #load ref data and build single cell experiments
    setwd(path_ref_data)
    
    ref_file_symbol <- load_data(path_ref_data, ref_file_symbol_name,'ref_symbol',opt_test,opt_load_from_loom)
    ref_map <- load_data(path_ref_data, ref_file_map_name, 'map',opt_test,opt_load_from_loom)
    
    sce_ref <- SingleCellExperiment(assays=list(normcounts=as.matrix(ref_file_symbol)), colData = ref_map)
    logcounts(sce_ref) <- log2(normcounts(sce_ref) + 1)
    # use gene names as feature symbols
    rowData(sce_ref)$feature_symbol <- rownames(sce_ref)
    # remove features with duplicated names
    sce_ref <- sce_ref[!duplicated(rownames(sce_ref)), ]
    mem_used()
    
    if (m=='highest_dropout_rate' & opt_determine_HDR_genes==TRUE){
      genes_selected <- get_genes_with_highest_dropout_rate(sce_ref,path_results_figures,ng)
    } else{
      #rowData(sce_ref)$scmap_features <- ... #set all to true
      rowData(sce_ref)$scmap_features <- logical(length=dim(sce_ref)[1])==FALSE 
    }
    
    query_file_name <- query_filelist[as.numeric(args[1])]
    print(query_file_name)
    
	  print(path_query_data)
    setwd(path_query_data)
    query_data <- load_data(path_query_data,query_file_name,'query',opt_test,opt_load_from_loom)
    if (opt_test==TRUE){
      query_GT <- load_data(path_query_data,query_file_GT,'map',opt_test,opt_load_from_loom)
      sce_query <- SingleCellExperiment(assays=list(normcounts=as.matrix(query_data)), colData = query_GT)
    } else{
      sce_query <- SingleCellExperiment(assays=list(normcounts=as.matrix(query_data)))
    }
    logcounts(sce_query) <- log2(normcounts(sce_query) + 1)
    rowData(sce_query)$feature_symbol <- rownames(sce_query)
    sce_query <- sce_query[!duplicated(rownames(sce_query)), ]
    
    if (opt_test==TRUE){
      add_str_query_file_name <- ""
    } else{ 
      add_str_query_file_name <- paste("_",paste(substr(query_file_name,1,3),substr(query_file_name,nchar(query_file_name)-5,nchar(query_file_name)-4),sep=""),sep="")
    }
    
    add_str <- paste(paste(paste(ng,sep=""),'_',sep=""),m,sep="")
    
    for (o in opt_scmap){
      if (o=='cluster'){
        #scmap-cluster
        sce_ref <- indexCluster(sce_ref)
        
        #some plots:
        png(file = paste(paste(paste(path_results_figures,'heatmap_cluster_',sep=""),add_str,sep=""),'.png',sep=""),width=800, height=1000)
        pheatmap(as.matrix(metadata(sce_ref)$scmap_cluster_index), show_rownames = FALSE)
        dev.off()
        
        scmapCluster_results <- scmapCluster(projection = sce_query, index_list = list(ABM = metadata(sce_ref)$scmap_cluster_index),threshold =0.5)
        if (opt_test==TRUE){
          png(file = paste(paste(paste(path_results_figures,'sankey_cluster_',sep=""),add_str,sep=""),'.png',sep=""),width=350, height=600)
          plot(getSankey(colData(sce_query)$cell_type1, scmapCluster_results$scmap_cluster_labs[,"ABM"],plot_width = 1200, plot_height = 600))
          dev.off()
        }
        
        R_cluster_ann_CT <- scmapCluster_results$scmap_cluster_labs # for all cells in test data set CT label
        R_cluster_score <- scmapCluster_results$scmap_cluster_siml#similarity score
        #head(scmapCluster_results$combined_labs) # if several ref data sets
        
        #TO DO: debug, write as function and call for other 2 methods too
        save_results_from_scmap(query_file_name, R_cluster_score, R_cluster_ann_CT, n_CT, o,ref_data_folder, opt_save_as_loom, add_str_query_file_name,opt_doublet_version, add_str)
          
        if (opt_test==TRUE){
          settings <- c(settings, add_str)
          kappa_val <- c(kappa_val,kappa2(data.frame(colData(sce_query)$cell_type1,scmapCluster_results$scmap_cluster_labs[,"ABM"])[scmapCluster_results$scmap_cluster_labs[,"ABM"]!="unassigned",])$value)
          p_unassigned <- c(p_unassigned,length(which(scmapCluster_results$scmap_cluster_labs[,"ABM"]=="unassigned"))/ncol(sce_query)*100)
        } 
      } else {
        #o=='cell':
        #scmap-cell
        
        #0. reduce ref data to only the genes present in the query data
        list_common_genes <- intersect(rowData(sce_ref)$feature_symbol, rowData(sce_query)$feature_symbol)
        sce_ref_red <- sce_ref[rowData(sce_ref)$feature_symbol %in% list_common_genes, ]
        sce_query_red <- sce_query[rowData(sce_query)$feature_symbol %in% list_common_genes, ]
        
        #1. Create index, prepare ref data
        set.seed(128746) # for stochasticity in k-means step
        sce_ref_red <- indexCell(sce_ref_red)
        
        #2. projection of test data
        scmapCell_results <- scmapCell(projection = sce_query_red, list(ABM = metadata(sce_ref_red)$scmap_cell_index))
        R_cell_ann_cellID <- scmapCell_results$ABM$cells
        R_cell_score <- scmapCell_results$ABM$similarities
        
        save_results_from_scmap(query_file_name, R_cell_score, R_cell_ann_cellID, n_CT, o,ref_data_folder, opt_save_as_loom, add_str_query_file_name,opt_doublet_version, add_str)
        
        if (o=='cell2cluster'){
          #scmap-cell2cluster
          # annotates the cells from the projection dataset using the labels of the reference --> looks at top 3 nearest neighbors and if they all belong
          # to the same cluster in the ref and their max similarity is higher than 0.5 --> projection cell is assigned to the corresponding reference cluster
          scmapCell_clusters <-scmapCell2Cluster(scmapCell_results,list(colData(sce_ref_red)$cell_type1))
          R_cell2cluster_ann_CT <- scmapCell_clusters$scmap_cluster_labs
          R_cell2cluster_score <- scmapCell_clusters$scmap_cluster_siml
          
          save_results_from_scmap(query_file_name, R_cell2cluster_score, R_cell2cluster_ann_CT, n_CT, o,ref_data_folder, opt_save_as_loom, add_str_query_file_name,opt_doublet_version, add_str)
        }
      }
    }
  }
}

stopCluster(cl)

if (opt_test==TRUE){
  Quality_cluster <- data.frame(settings, kappa_val,p_unassigned)
  write.csv(Quality_cluster,paste(paste(paste(paste(path_query_data,"matrix_",sep=""),n_cells_query,sep=""),"_cells/scmap/",sep=""),"Quality_cluster",sep=""))
}
