#author: Lisa Bast
#date: 18.01.21
#version: 0.0.2
##helper functions for scmap

#determine top k genes with highest dropout rate and store these in a csv file
get_genes_with_highest_dropout_rate <- function(sce_ref,path_results_figures,ng){
  png(file=paste(paste(path_results_figures,'selection_',ng,sep=""),'_features_based_on_dropout.png',sep=""), width=600, height=350)
  sce_ref <- selectFeatures(sce_ref, n_features = ng, suppress_plot = FALSE)
  dev.off()
  
  #table(rowData(sce_ref)$scmap_features)
  genes_selected<-rownames(sce_ref)[rowData(sce_ref)$scmap_features]
  write.csv(genes_selected,paste(paste(paste(path_ref_data,"list_genes_selected_high_dropout_",sep=""),ng,sep=""),".csv",sep=""))
  return(genes_selected)
}

#load data for cell type annotation with scmap
load_data <- function(path, filename, file_content, opt_test, opt_loom) {
  print(path)
  setwd(path)
  #if (file_content=='ref_symbol'){
  #  D <- read.csv.raw(file = filename, nrows=2000)
  #} else {
  #  D <- read.csv.raw(file = filename)
  #}
  if ((opt_loom==TRUE)&(file_content=='query')){
    D_loom <- connect(filename,mode="r",skip.validate =TRUE)
    D <- t(D_loom[["matrix"]][, ])
    row_names <- D_loom[["row_attrs/Gene"]][]
    #close loom file connection:
    D_loom$close_all()
  }
  else{
    D <- read.csv.raw(file = filename)
  }
  if (file_content=='query'){
    if (opt_loom==FALSE){
      row_names <- D[,1]
      c_names <- colnames(D)
      D <- D[c_names[c_names!=""]]
    } 
    row.names(D) <- row_names
    # remove rows of gene duplicates
    if (opt_test==FALSE){
      D <- D[which(duplicated(row_names)==FALSE),]
    }
    if ((opt_test==FALSE)&(opt_loom==TRUE)){
      D<-D*10/10
    }
  } else if (file_content=='ref_symbol'){
    row.names(D) <- D$symbol
    c_names <- colnames(D)
    D <-D[c_names[c_names!="symbol" & c_names!=""]]
  } else if (file_content=='map'){
    row.names(D) <- D$sample
    D <- D[c("sample","cell type")]
    colnames(D)<-c(" ","cell_type1")
    D <-D["cell_type1"]
  } else{
    stop('The specified file content type does not exist!')
  }
  return(D)
}

#save result from cell type annotation with scmap
save_results_from_scmap <- function(filename, R_cluster_score, R_cluster_ann_CT, n_CT, o,ref_data_folder, opt_save_as_loom, add_str_query_file_name,opt_doublet_version, add_str) {
  if (opt_save_as_loom==TRUE){
    #connect to loom with read and write access
    lfile <- connect(filename, mode = "r+", skip.validate=TRUE)
    #directly integrate R_cluster_ann_CT into loom file metadata
    #directly integrate R_cluster_score into loom file metadata
    new_ca_name_score = paste0(paste0(paste0(paste0(paste0(paste0('CT_ann_score_',str_sub(ref_data_folder,1,nchar(ref_data_folder)-1)),'_scmap_'),o),'_'),as.character(n_CT)),'_CTs')       
    new_ca_name_ann = paste0(paste0(paste0(paste0(paste0(paste0('CT_ann_',str_sub(ref_data_folder,1,nchar(ref_data_folder)-1)),'_scmap_'),o),'_'),as.character(n_CT)),'_CTs') 
    cat(new_ca_name_score)
    cat(new_ca_name_ann)
    #replace na with string na
    R_cluster_ann_CT[,'ABM'][is.na(R_cluster_ann_CT[,'ABM'])]='nan'
    col_attrs <- list(R_cluster_score[,'ABM'],R_cluster_ann_CT[,'ABM'])
    names(col_attrs) <- c(new_ca_name_score,new_ca_name_ann)
    lfile$add.col.attribute(col_attrs,overwrite=TRUE)
    #close loom file
    lfile$close_all()
  } 
  else{
    filenames <- c(paste0(paste0(paste0(paste0(paste0(paste0("R_",o),"_ann_CT_"),add_str),add_str_query_file_name),opt_doublet_version),".csv"),paste0(paste0(paste0(paste0(paste0(paste0("R_",o),"_score_"),add_str),add_str_query_file_name),opt_doublet_version),".csv"))
    
    write.csv(R_cluster_ann_CT,paste(path_results_files,filenames[1],sep=""))
    write.csv(R_cluster_score,paste(path_results_files,filenames[2],sep=""))
  }
}