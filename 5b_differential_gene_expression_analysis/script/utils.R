#helper functions to load data, put together and store results, plot or load intermediate results
#author: Lisa Bast
#date: 27.04.2021
#version: 0.0.3


#settings:
get_paths <- function(main_path,main_project_path,data_folder){
  paths <- c()
  paths$sample_info_data <- paste0(main_project_path,"/2_alignment/output/")
  paths$QC_data <- paste0(main_project_path,"/3_quality_control/output/")
  paths$data <- paste0(main_path,"/output/",data_folder)
  paths$results <- paste0(main_path,"/output/DEGs/")
  paths$results_RUV <- paste0(main_path,'/output/RUV/')
  paths$xy_chr_list <- paste0(main_project_path,'/4_data_integration_and_cell_type_annotation/output/')
  paths$CT_spec_file<-paste0(main_path,'/data/')
  #make sure path exists
  dir.create(paths$results, recursive=TRUE, showWarnings = FALSE)
  return(paths)
}

get_settings <- function(opt_input_data){
  settings <- c()
  if (opt_input_data=="GRN_modules_gprofiler_result"){
    settings$opt_up_and_down_sep = FALSE #
  } else if (opt_input_data=="Patient_clustering_genes_GSA_result" | opt_input_data=="GRN_modules_GSA_result" | opt_input_data=="all_GRN_modules_GSA_result"){
    settings$opt_up_and_down_sep = FALSE #
  } else{
    settings$opt_up_and_down_sep = TRUE#FALSE #
  }
  
  opt_with_MT_genes <- TRUE#FALSE # with or without mito dna genes
  opt_new_profiler_settings <- FALSE #TRUE
  settings$opt_metacell_settings <- "50_13_250" #"30_6_120"#
  
  if (opt_input_data == "Patient_clustering_genes_GSA_result"){
    settings$n_clusters_range = c('SCZ_cases')#,'whole_dataset','controls')
    settings$min_overlap_range <- c(3)#c(3,5,10)
    settings$background_folders <- c('cortex/','SnRNAseqDet/')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$CT_groups <- c('all')
  } else if (opt_input_data =="Patient_clustering_genes_GSEA_result"){
    settings$n_clusters_range = c('controls','SCZ_cases','whole_dataset')
    settings$min_overlap_range <- c(0)
    settings$background_folders <- c('none/')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$CT_groups <- c('all')
  } else if (opt_input_data == "DEG_GSA_result"){ 
    settings$n_clusters_range <- c(16)#,37,3)
    settings$sub_folders <- c('alpha3/')
    settings$min_overlap_range <- c(3)#c(3,5,10)
    settings$background_folders <- c('cortex/')#,'cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$CT_groups <- c("Excitatory_Layer_5_6_IT_neurons_I","Excitatory_Layer_5_6_IT_neurons_II",
                   "Excitatory_Layer_5_6_CT_and_NP_neurons","Excitatory_Layer_2_3_IT_neurons_I",
                   "Excitatory_Layer_2_3_IT_neurons_II", "Excitatory_Layer_3_4_IT_neurons","Excitatory_Layer_3_6_IT_neurons",
                   "Astrocytes","Endothelial_and_mural_cells","Microglial_cells","Oligodendrocyte_progenitor_cells","Oligodendrocytes",
                   "Inhibitory_LAMP5_neurons","Inhibitory_PVALB_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons")                     
    #old definition:                                        
    #CT_groups <- c('all','Astrocytes', "Microglial_cells","Oligodendrocyte_progenitor_cells",'Neurons')
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){ 
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    settings$sub_folders <- c('norm_basic/')
    settings$min_overlap_range <- c(0)
    settings$background_folders <- c('none/')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$opt_modules_separately <- TRUE
    settings$CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
    settings$n_clusters_range <- c(15)
    settings$sub_folders <- c('alpha3/')
    settings$min_overlap_range <- c(3)#c(3,5,10)
    settings$background_folders <- c('cortex/')#,'cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    settings$opt_sim_based_on <- 'size'
    settings$CT_groups <- c('all')
  } else if(opt_input_data == "all_GRN_modules_gprofiler_result"){
    settings$n_clusters_range <- c(15)
    settings$sub_folders <- c('alpha3/')
    settings$min_overlap_range <- c(3)#c(3,5,10)
    settings$background_folders <- c('cortex/')#,'cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    settings$opt_sim_based_on <- 'size'
    settings$CT_groups <- c('all')
  } else if (opt_input_data == "all_GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    settings$sub_folders <- c('GSA/')
    settings$min_overlap_range <- c('')
    settings$background_folders <- c('')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$opt_modules_separately <- TRUE
    settings$CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if (opt_input_data == "GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    settings$sub_folders <- c('GSA/')
    settings$min_overlap_range <- c('')
    settings$background_folders <- c('')
    settings$opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    settings$opt_modules_separately <- FALSE
    settings$CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data =="longread_results"){
    settings$n_clusters_range <- c('dummy')
    settings$sub_folders <- c('')
    settings$min_overlap_range <- c('')
    settings$background_folders <- c('')
    settings$opt_sim_based_on <- 'size'
    settings$CT_groups <- c('all')
  }
  return(settings)
}

get_metadata <- function(data_all, paths, filename_data_sample_info, filename_data_ancestry, filename_QC_data, opt_load_metadata, opt_scale_and_center_metadata,opt_ancestry_binary){

  if (opt_scale_and_center_metadata==TRUE){
    add_str_scaled_and_centered <- "_scaled_and_centered"
  } else{
    add_str_scaled_and_centered <- ""
  }
  if (opt_load_metadata==FALSE){
    #vars to add: donor_ID_universal, donor_ID_internal_6, donor_ID_internal_8, pH_of_CSF, source 
    data_sample_info <- read_excel(paste0(paths$sample_info_data,filename_data_sample_info), sheet = "basic_donor_information",col_names = TRUE)[,c("donor_ID_universal","donor_ID_internal_8","genotype_IID","source","gender","age")]#3,7,9,10,11)]
    

    bool_keep <- data_sample_info$donor_ID_internal_8 %in% data_all$colData$donor_ID_python
    data_sample_info <- data_sample_info[bool_keep, ]
    names(data_sample_info)[names(data_sample_info)=="donor_ID_internal_8"] <- "donor_ID_python"
    names(data_sample_info)[names(data_sample_info)=="genotype_IID"] <- "donor_ID_internal_6"
    names(data_sample_info)[names(data_sample_info)=="gender"] <- "Sex"
    names(data_sample_info)[names(data_sample_info)=="age"] <- "Age"
    #read file for ancestry (binary)
    data_ancestry <- read.table(file = paste0(paths$sample_info_data,filename_data_ancestry), sep = '\t', header = TRUE)[,c("IID", "genotype_if_EUR", "genotype_batch")]
    names(data_ancestry)[names(data_ancestry)=="genotype_if_EUR"] <- "ancestry_european"
    #split data set according to IDs
    data_ancestry_ID_uni <- data_ancestry[data_ancestry$genotype_batch=="brn2",][,c(1,2)]
    data_ancestry_ID_D6 <- data_ancestry[data_ancestry$genotype_batch=="brn1",][,c(1,2)]
    #rename columns:
    names(data_ancestry_ID_uni)[names(data_ancestry_ID_uni)=="IID"] <- "donor_ID_universal"
    names(data_ancestry_ID_D6)[names(data_ancestry_ID_D6)=="IID"] <- "donor_ID_internal_6"

    
    #merge:
    data_sample_info_merged_1 <- merge(data_sample_info,data_ancestry_ID_uni,by="donor_ID_universal")
    data_sample_info_merged_2 <- merge(data_sample_info,data_ancestry_ID_D6,by="donor_ID_internal_6")
    data_sample_info_merged <- rbind(data_sample_info_merged_1,data_sample_info_merged_2)
    #make sure all donors are in df, even if ancestry is missing:
    data_sample_info_merged <- merge(data_sample_info_merged,data_sample_info,all.y=T)
    
    #data_sample_info_merged <- merge(data_sample_info,data_sample_info_merged,by=c("donor_ID_universal","donor_ID_internal_6","donor_ID_python","source","pH_of_CSF"),all.x=T,all=T)
    #remove IDs not used:
    data_sample_info_merged <- subset(data_sample_info_merged, select=-c(donor_ID_universal,donor_ID_internal_6))

    data_QC <- read_excel(paste0(paths$QC_data,filename_QC_data))[,c(3:20)]#21#contains 83 samples, no NAs, S78, S66 missing
    data_QC <- subset(data_QC,select=-c(mean_reads_per_umi,Group,Library,age,p_cells_del_filt))
    
    names(data_QC)[names(data_QC)=="gender"] <- "Sex"
    names(data_QC)[names(data_QC)=="sex"] <- "Sex"
    #names(data_QC)[names(data_QC)=="age"] <- "Age"
    
    #vars_to_add <- c("p_valid_barcodes","p_sequencing_saturation","p_genome_not_gene","p_mapped_reads","p_unmapped_reads","mean_counts_per_barcode","std_counts_per_barcode","median_gpc")  
    colData <- data_all$colData[,c("donor_ID_python", "p_cells_del_filt", "mean_reads_per_umi", "Group", "Age", "Library")]
    
    
    metadata_merged <- merge(colData, data_QC[data_QC$donor_ID_python %in% colData$donor_ID_python, ], by=c("donor_ID_python"), all=T)
    metadata_merged <- subset(metadata_merged,select=-c(Sex,Age))

    metadata_merged_2 <- merge(data_sample_info_merged[data_sample_info_merged$donor_ID_python %in% colData$donor_ID_python, ], metadata_merged[metadata_merged$donor_ID_python %in% colData$donor_ID_python, ], by=c("donor_ID_python"), all=T)  

    #reorder metadata according to how samples appear in data_all object:
    metadata_sorted <- metadata_merged_2[match(data_all$colData$donor_ID_python, metadata_merged_2$donor_ID_python), ]  
    
    ## update data_all$colData (integrate metadata_sorted)
    #data_all$colData <- metadata_sorted
    
    #which variables in metadata_sorted contain missing values?
    NA_idx <- which(apply(metadata_sorted, 2, function(x) any(is.na(x)))==TRUE)
    num_NAs <- apply(metadata_sorted, 2, function(x) sum(is.na(x)==TRUE))
    names(NA_idx)<-NULL
    
    if (paths$results==''){
      png(file=paste0(paths$results,"missing_data_in_metadata_per_variable.png"),width=600, height=550, units='mm', res=300)
      par(mfrow=c(1,1),mar=c(1,1,1,1)) 
      mice_plot <- aggr(metadata_sorted, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(metadata_sorted), cex.axis=0.45, gap=1, ylab=c("Missing data","Pattern"))
      dev.off()
      graphics.off()
    }
    #impute missing values
    #imputation method: classification and regression trees (stochastic method)
    metadata_imputed_tmp <- mice(metadata_sorted, m=5, maxit = 500, method = 'cart', seed = 500)
    metadata_imputed <- complete(metadata_imputed_tmp,1)

    #make sure metadata columns are either factors or numeric
    metadata_imputed$source <- as.factor(metadata_imputed$source)
    if (opt_ancestry_binary == T){
      metadata_imputed$ancestry_european <- as.factor(metadata_imputed$ancestry_european)
    }
    metadata_imputed$Group <- as.factor(metadata_imputed$Group)
    
    metadata_imputed$Library <- as.factor(metadata_imputed$Library)
    metadata_imputed$Sex <- as.factor(metadata_imputed$Sex)
    
    #if necessary scale and center metadata
    if (opt_scale_and_center_metadata==TRUE){
      cols_numeric = c("Age","mean_reads_per_umi", "p_cells_del_filt","num_umi","num_reads","p_valid_barcodes","p_sequencing_saturation","p_genome_not_gene",
                       "p_mapped_reads","p_unmapped_reads","PMI_h","mean_counts_per_barcode","std_counts_per_barcode","median_gpc")
      bool_col_to_be_trans <- rep(FALSE,length(colnames(metadata_imputed)))
      bool_col_to_be_trans[colnames(metadata_imputed) %in% cols_numeric] <- TRUE
      metadata_imputed_scaled <- metadata_imputed
      for (i in seq(1,length(bool_col_to_be_trans))){
        if (bool_col_to_be_trans[i]==TRUE){
          metadata_imputed_scaled[,i]<-(metadata_imputed[,i]-mean(metadata_imputed[,i]))/sqrt(var(metadata_imputed[,i]))
        }
      }
      #plot PCA and loadings
      if (paths$results!=''){
        plot_PCA_of_metadata(metadata_imputed_scaled,cols_numeric,paths$results)
      }
      metadata_imputed <- metadata_imputed_scaled
    }
    
    #save as file:
    write.csv(metadata_imputed,paste0(paths$sample_info_data,'metadata_imputed_cellranger',add_str_scaled_and_centered,'.csv'),row.names = FALSE)
  } else{
    #make sure the scaled and centered data is read in case this was specified:
    metadata_imputed_complete <- read.csv(paste0(paths$sample_info_data,'metadata_imputed_cellranger',add_str_scaled_and_centered,'.csv'))
    # make sure the right samples are in metadata (same as in data$CC)
    metadata_imputed <- metadata_imputed_complete[metadata_imputed_complete$donor_ID_python %in% data_all$colData$donor_ID_python, ]
  }
  
  return(metadata_imputed)
}

get_cell_Type_name_from_filename <- function(countData_filenames,file_id,opt_aggregation){

  file_name_pieces <- str_split(countData_filenames[file_id],'_')
  celltype_aggr <- str_split(countData_filenames[file_id],paste0('_',opt_aggregation))[[1]][1]
  celltype <- str_split(celltype_aggr[1],'_filtered_')[[1]][2]
  
  return(celltype)
}

#function to load aggregated data
load_count_data <- function(countData_filenames,colData_filenames,file_id,opt_aggregation,paths){
  

  setwd(paths$data)
  celltype <- get_cell_Type_name_from_filename(countData_filenames,file_id,opt_aggregation)
  file_name_pieces <- str_split(countData_filenames[file_id],'_')
  results_subfolder_raw <- paste0(do.call(paste,c(as.list(file_name_pieces[[1]][9:length(file_name_pieces[[1]])]),sep="_")))
  results_subfolder <- paste0(substr(results_subfolder_raw,1,nchar(results_subfolder_raw)-4),'/')
  
  #countData: matrix of non-negative integers
  CC<-read.csv.raw(file = countData_filenames[file_id])
  #colData: dataframe with at least a sinlge column; rows of colData correspond to columns of countData
  colData<-data.frame(read.csv.raw(file = colData_filenames[file_id]))
  colData$num_umi <- scale(colData$num_umi, center = TRUE, scale = TRUE)
  
  
  #make gene names unique or use unique Ensembl gene IDs:
  if ("Accession" %in% colnames(CC)){
    row.names(CC) = paste(CC[,which(colnames(CC)=="Gene")],CC[,which(colnames(CC)=="Accession")],sep='_')
    col_sel_start <- which(colnames(CC)=="Accession")+1
    CC <- CC[,col_sel_start:dim(CC)[2]]
  } else {
    gene_list <- CC[,which(colnames(CC)=="Gene")]
    duplicated_genes <- gene_list[duplicated(gene_list)]
    for (dg_id in seq(1,length(duplicated_genes))){
      CC[CC["Gene"]==duplicated_genes[dg_id],1] = paste(CC[CC["Gene"]==duplicated_genes[dg_id],1],c("_1","_2"),sep='')
    }
    row.names(CC) <- CC[,which(colnames(CC)=="Gene")]
    col_sel_start <- which(colnames(CC)=="Gene")+1
    CC <- CC[,col_sel_start:dim(CC)[2]]
  }
  
  #row.names(colData) <- colnames(CC)
  colnames(CC) <-NULL
  
  #remove X and Y chromosome genes
  xy_chr_genes <- read.csv.raw(paste0(paths$xy_chr_list,"XYchr_genes.csv"))
  bool_xy_chr <- rownames(CC) %in% paste(xy_chr_genes$Gene,xy_chr_genes$Accession,sep='_')
  CC <- CC[bool_xy_chr==FALSE, ]
  
  return(list("data_obj", colData=colData, CC=CC, celltype=celltype, results_subfolder=results_subfolder))
}

#randomize order of group labels:
get_randomized_group_labels <- function(group_labels_ordered, path_results, i){
  file_name <- paste0(paste0(paste0(path_results,"group_labels_random_"),as.character(i)),".csv")
  #check if this randomization already exists:
  if (file.exists(file_name)){
    group_labels_random_df <- read.csv(file = file_name)
    # make sure the vector has the same length as group_labels_ordered:
    group_labels_random_full <- group_labels_random_df$group_labels_random 
    group_labels_random <- group_labels_random_full[1:length(group_labels_ordered)]
  } else {
    group_labels_random <- sample(group_labels_ordered)
    #save the current random group labels:
    write.csv(data.frame(group_labels_random,group_labels_ordered), file=file_name, row.names = FALSE)
  }
  return(group_labels_random)
}

#filter results dataframe for significant genes and add settings as columns
get_significant_results <- function(R,res,alpha_val,sfe_type,shrinkage_type,celltype){
  #remove nans:
  if (length(res$baseMean[!is.na(res$padj)])>0){
    if (shrinkage_type=='none'){
      res_nan_free <- data.frame(baseMean=res$baseMean[!is.na(res$padj)],
                                 log2FoldChange=res$log2FoldChange[!is.na(res$padj)],
                                 lfcSE=res$lfcSE[!is.na(res$padj)],
                                 stat=res$stat[!is.na(res$padj)],
                                 pvalue=res$pvalue[!is.na(res$padj)],
                                 padj=res$padj[!is.na(res$padj)],stringsAsFactors=FALSE)
    } else{
      res_nan_free <- data.frame(baseMean=res$baseMean[!is.na(res$padj)],
                                 log2FoldChange=res$log2FoldChange[!is.na(res$padj)],
                                 lfcSE=res$lfcSE[!is.na(res$padj)],
                                 stat=rep(NaN,times=length(res$lfcSE[!is.na(res$padj)])),
                                 pvalue=res$pvalue[!is.na(res$padj)],
                                 padj=res$padj[!is.na(res$padj)],stringsAsFactors=FALSE)
    }
    rownames(res_nan_free) <- rownames(res)[!is.na(res$padj)]
    if (length(res_nan_free$baseMean[res_nan_free$padj<=alpha_val])>0){
      new_rows <- data.frame(baseMean=res_nan_free$baseMean[res_nan_free$padj<=alpha_val],
                             log2FoldChange=res_nan_free$log2FoldChange[res_nan_free$padj<=alpha_val],
                             lfcSE=res_nan_free$lfcSE[res_nan_free$padj<=alpha_val],
                             stat=res_nan_free$stat[res_nan_free$padj<=alpha_val],
                             pvalue=res_nan_free$pvalue[res_nan_free$padj<=alpha_val],
                             padj=res_nan_free$padj[res_nan_free$padj<=alpha_val],
                             size_factor_estimation=rep(sfe_type,sum(res_nan_free$padj<=alpha_val, na.rm=TRUE)),
                             shrinkage=rep(shrinkage_type,sum(res_nan_free$padj<=alpha_val, na.rm=TRUE)),
                             celltype=rep(celltype,sum(res_nan_free$padj<=alpha_val, na.rm=TRUE)),stringsAsFactors=FALSE)
      rownames(new_rows)<-rownames(res_nan_free)[res_nan_free$padj<=alpha_val]
      R <-rbind(R,new_rows)
    }
  }
  return(R)
}

plot_RUV_factor_correlation_with_metadata <- function(RUV_factors,metadata,results_path){
  
  RUV_factor_names <- colnames(RUV_factors)
  corr_types <- c("pearson","spearman")
  #merge data frames
  RUV_factors$donor_ID_python <- rownames(RUV_factors)
  Data <- merge(RUV_factors,metadata,by = "donor_ID_python")
  #continuous variables:
  var_cont <- c("Age" ,"PMI_h","p_cells_del_filt","mean_reads_per_umi","num_umi","num_reads","p_valid_barcodes","p_sequencing_saturation","p_genome_not_gene", "p_mapped_reads","p_unmapped_reads","mean_counts_per_barcode","std_counts_per_barcode","median_gpc")
  n <- length(RUV_factor_names)
  m <- length(var_cont)
  for (corr_method in corr_types){
    M_cor <- matrix(0,n,m)
    for (f_id in seq(1,n)){
      for (m_id in seq(1,m)){
        print(var_cont[m_id])
        M_cor[f_id,m_id] = cor(Data[RUV_factor_names[f_id]],Data[var_cont[m_id]],method=corr_method)
      }
    }
    rownames(M_cor) <- RUV_factor_names
    colnames(M_cor) <- var_cont
    #plot corr as ellipses
    pdf(file = paste0(results_path,as.character(length(RUV_factor_names)),"_RUV_factor_",corr_method,"_correlation_wih_metadata_continuous.pdf")) # The height of the plot in inches
    print(corrplot(M_cor,method="ellipse",addCoef.col = 'black',tl.col = 'black'))
    dev.off()
  }
  
  #categorical variables:
  browser()
  var_cat <- c("Group","source","Sex","Library") #"ancestry_european"
  m <- length(var_cat)
  p <- list()
  for (f_id in seq(1,n)){
    eval(parse(text=paste0("Data %>% pivot_longer(cols = var_cat) %>% 
      ggplot(aes(x=factor(value),y=",RUV_factor_names[f_id],")) + geom_boxplot() + facet_wrap(.~name,scales='free_x') + theme_bw() + xlab('metadata variables') + geom_jitter(pch=21,position = position_jitter(width = 0.2),show.legend = FALSE,alpha = 0.5)
    ggsave(paste0(results_path,RUV_factor_names[f_id],'_boxplots_for_categorical_metadata.pdf'))")))
  }
  #pdf(file = paste0(results_path,as.character(length(RUV_factor_names)),"_RUV_factor_boxplots_for_categorical_metadata.pdf"))
  #do.call(grid.arrange,p)
  #dev.off()
}

#remove genes with low counts and store data in DDS object:
get_filtered_data <- function(CC,colData,design_formula,path_xy_chr_list,p_samples_larger_10_counts){
  #remove X and Y chromosome genes
  setwd(path_xy_chr_list)
  xy_chr_genes <- read.csv.raw("XYchr_genes.csv")
  bool_xy_chr <- rownames(CC) %in% paste(xy_chr_genes$Gene,xy_chr_genes$Accession,sep='_')
  CC <- CC[bool_xy_chr==FALSE, ]

  DDS <- DESeqDataSetFromMatrix(countData=round(CC),colData,formula(eval(parse(text=design_formula))),tidy=FALSE,ignoreRank=FALSE)
  if (opt_filter_genes_with_low_counts==TRUE){
    # pre-filtering: remove rows in which there are very few reads:
    keep <- rowSums(counts(DDS)) >= 10
    DDS <- DDS[keep,]
    # remove low count genes:
    #at least p% of samples need to have more than 10 counts
    n_samples_min <- dim(DDS)[2]*(p_samples_larger_10_counts/100)
    DDS <- DDS[rowSums(counts(DDS)>=10)>n_samples_min,]
  }
  return(DDS)
}

get_indices_for_random_dataset_halves <- function(DDS_RUV,min_p_donors_groups){
  n_samples <- dim(counts(DDS_RUV))[2]
  first <- sort(sample(1:n_samples, floor(n_samples/2), replace=F))
  # in case less than 13 individuals in one of the groups draw random numbers again
  n_SCZ <- sum(DDS_RUV$Group=='SCZ')
  n_CTRL <- sum(DDS_RUV$Group=='CTRL')
  while (sum(DDS_RUV$Group[first]=='SCZ')<round(min_p_donors_groups*n_SCZ/2) | sum(DDS_RUV$Group[first]=='CTRL')<round(min_p_donors_groups*n_CTRL/2)){ 
    first <- sort(sample(1:n_samples, floor(n_samples/2), replace=F))
  }
  second <-which((seq(1,n_samples) %in% first) == FALSE)
  idx_rnd <- list("first"=first, "second"=second)
  # data_RUV$CC_1st_half <- data_RUV$CC[,idx_rnd_1st]
  # data_RUV$CC_2nd_half <- data_RUV$CC[,idx_rnd_2nd]
  # data_RUV$metadata_1st_half <- data_RUV$metadata[idx_rnd_1st,]
  # data_RUV$metadata_2nd_half <- data_RUV$metadata[idx_rnd_2nd,]
  return(idx_rnd)
}


#decide which sfe type based on data
#determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
get_sfe_type <- function(DDS){
  n_genes_pos_in_all_samples <- sum(rowSums(counts(DDS)>0)==dim(DDS)[2])
  #if >1000 genes --> use poscounts, otherwise use ratio
  if (n_genes_pos_in_all_samples>1000){
    sfe_type <- c("poscounts")
  } else{
    sfe_type <- c("ratio")
  }
  return(sfe_type)
}

#get and create results folder
get_results_folder<-function(path_res,design_interactions_str,results_subfolder,i,opt_randomize_group_labels){
  

  if (opt_randomize_group_labels==TRUE){
    path_res <- paste0(path_res,"randomized_group_labels_",as.character(i),"/")
  }
  path_res <- paste0(path_res,design_interactions_str,"/")
  #dir.create(path_res_R,recursive=TRUE, showWarnings = FALSE)
  if (results_subfolder!='all'){
    path_res <- paste0(path_res,results_subfolder,'/')
  } 
  
  dir.create(path_res,recursive=TRUE, showWarnings = FALSE)
  return(path_res)
}

#sum aggregated data across cell types
#date: 13.10.2021
load_count_data_aggregated_across_celltypes <- function(paths, path_code, opt_aggregation, add_str_pagoda, opt_downsampled, i){

  setwd(paths$data)
  if (opt_downsampled==TRUE){
    file_str_end <- paste0("_ds_",i,".csv")
    countData_search_str <- paste0("Aggr_counts",add_str_pagoda,"_TH_and_D_adj_filtered_all_",opt_aggregation,file_str_end)
    colData_search_str <- paste0("G_info",add_str_pagoda,"_TH_and_D_adj_filtered_all_",opt_aggregation,file_str_end)
  } else{
    file_str_end <- "*.csv"
    countData_search_str <- paste0("Agg_counts_",add_str_pagoda,"_TH_and_D_adj_filtered_all_",opt_aggregation,file_str_end)
    colData_search_str <- paste0("G_info_",add_str_pagoda,"_TH_and_D_adj_filtered_all_",opt_aggregation,file_str_end)
  }
  
  print(countData_search_str)
  print(colData_search_str)
  
  countData_filenames <- Sys.glob(countData_search_str)
  colData_filenames <- Sys.glob(colData_search_str)
  
  setwd(path_code)
  data_all <- load_count_data(countData_filenames,colData_filenames,1,opt_aggregation,paths)

  #aggregate data across cell types:
  #for (file_id in seq(1,length(countData_filenames))){
  #  #load data:
  #  data <- load_count_data(countData_filenames,colData_filenames,file_id,opt_celltype_groups,opt_aggregation,paths)
  #  if (file_id==1){
  #    data_all <- data
  #    donors <- data_all$colData['donor_ID_python'][[1]]
  #  } else{
  #    data_all$CC[,1] <- data_all$CC[,1] + data$CC[,1]
  #  }
  #}
  
  #remove X and Y chromosome genes
  setwd(paths$xy_chr_list)
  xy_chr_genes <- read.csv.raw("XYchr_genes.csv")
  bool_xy_chr <- rownames(data_all$CC) %in% paste(xy_chr_genes$Gene,xy_chr_genes$Accession,sep='_')
  data_all$CC <- data_all$CC[bool_xy_chr==FALSE, ]
  
  #update metadata:
  data_all$celltype<-"all"
  data_all$results_subfolder<-"all"
  
  return(data_all)
}

#extract matrix of size factor normalized counts from DESeq2 data structure
get_sfe_normalized_data <- function(DDS){
  data_sfe_normalized <- counts(DDS)*t(replicate(dim(DDS)[1],DDS$sizeFactor))
}

#perform MARS and return variables left in MARS model
get_variables_selected_by_mars <- function(data_all,n_rep,path_results){
  
  n_genes <- dim(data_all$CC_sfe_normalized)[1]
  if (n_genes>1000){
    n_genes_target <- 1000
  } else if ((n_genes>=500) & (n_genes<1000)){
    n_genes_target <- 500
  } else{
    n_genes_target <- 1000
  }
  
  if (n_genes_target<=n_genes){
    for (i in seq(1,n_rep)){
      idx_rnd <- sample(1:n_genes, n_genes_target, replace=F)
      #check if always different genes are selected
      if (i==1){
        idx_rnd_1<-idx_rnd
      } else{
        if (identical(sort(idx_rnd),sort(idx_rnd_1))){
          print("gene set is identical with first one!")
          break
        }
      }
      
      Gene_set_i <- t(as.matrix(data_all$CC_sfe_normalized[idx_rnd,]))
      # scale gene expression:
      
      #rownames(Gene_set_i) <- NULL
      
      #metadata$Gene_i <- unlist(t(Gene_i), use.names=FALSE)
      #basic mars model:
      mars1 <- earth(
        Gene_set_i ~ .,
        data = subset(data_all$metadata, select = -c(donor_ID_python,Group)),
        trace=2
      );
      #show number of terms (vertical line)
      #somehow extract metadata variables that are in selection
      #plot(mars1, which=1)
      #summary(mars1) %>% .$coefficients
      if (i<=5){
        png(file=paste0(paste0(paste0(path_results,"Model_error_per_number_of_terms"),i),".png"),width=600, height=550, units='mm', res=300)
        par(mfrow=c(1,1),mar=c(1,1,1,1)) 
        plot(mars1, which = 1)
        dev.off()
        graphics.off()
      }
      
      selected_terms <- colnames(mars1$bx)[2:dim(mars1$bx)[2]]
      if ((i==1)&(i==1)){
        vars_selected <- matrix(FALSE,n_rep,length(rownames(mars1$modvars)))
        var_list<- list()
      }
      j <- 0
      for (var in rownames(mars1$modvars)){
        var_list[[var]]<- c()
        j=j+1
        for (term_sel in selected_terms){
          if (str_contains(term_sel,var)){
            print(term_sel)
            vars_selected[i,j] <- TRUE
            var_list[[var]]<- c(var_list[[var]],i)
            break
          }
        }
      }
      if (i==1){
        colnames(vars_selected)<-rownames(mars1$modvars)
      }
    }
  } else{
    vars_selected <- 303
  }
  return(vars_selected)
}

get_RUV_data <- function(DDS_RUV,res_RUV,opt_filter_high_CT_specificity_genes,path_CT_spec,filename_CT_spec,opt_data_halves,pval_TH,specificity_TH,min_p_donors_groups){
  if (opt_data_halves==TRUE){
    idx_rnd <- get_indices_for_random_dataset_halves(DDS_RUV,min_p_donors_groups)
  } else{
    idx_rnd <- c()
  }
  
  set_ini <- newSeqExpressionSet(counts(DDS_RUV))
  
  #only keep genes that are detected more than 5 times in at least 2 donors in both data halves:
  #idx  <- rowSums(counts(set_1st_half_ini) > 5) >= 2 | rowSums(counts(set_2nd_half_ini) > 5) >= 2
  idx  <- rowSums(counts(set_ini) > 5) >= 2 
  set_ini  <- set_ini[idx, ]
  
  set_ini <- betweenLaneNormalization(set_ini, which="upper")
  
  #determine subset of genes which should be used as negative controls in the estimation of factors of unwanted variation
  #genes far from being significant in initial deseq2 run:
  #not.sig <- rownames(res_RUV)[which(res_RUV$pvalue > pval_TH)]
  not.sig <- rownames(res_RUV)[which(res_RUV$pvalue > pval_TH)]
  
  empirical <- rownames(set_ini)[ rownames(set_ini) %in% not.sig ]
  
  print(paste0(length(empirical)," genes left in data set for RUV analysis before CT specificity filtering!"))
  
  #load cell type specificity scores for genes
  #filter out genes with high cell type specificity scores
  if (opt_filter_high_CT_specificity_genes==TRUE){
    load(paste0(path_CT_spec, filename_CT_spec))#object is called dat
    #sort genes by specificity
    dat_sorted <- dat[order(dat$specificity,decreasing=TRUE),]
    #create list of unique genes while preserving order
    genes_ordered_by_specificity <- unique(dat_sorted$ENSGID)
    #remove top (1-specificity_TH*100)% of genes
    genes_with_high_specificity <- genes_ordered_by_specificity[1:round(length(genes_ordered_by_specificity)*(1-specificity_TH))]
    
    #old implementation:
    #genes_with_high_specificity <- dat$ENSGID[dat$specificity>=quantile(dat$specificity,specificity_TH)]
    #does gene name in empirical end with ENSGID in dat? --> remove from empirical
    for (i in seq(1,length(genes_with_high_specificity))){
      empirical <- empirical[!endsWith(empirical,genes_with_high_specificity[i])]
    }
  }
  #to do: export list of genes used for RUV analyses
  #perform PCA
  
  print(paste0(length(empirical)," genes left in data set for RUV analysis!"))
  D <- list(idx_rnd=idx_rnd,set_ini=set_ini,empirical=empirical)
  return(D)
}

##for a particular data split: calculate log2FC of DESeq2 result between data halves for a range of RUV values
get_correlation_LFC_data_halves <- function(DDS, n_RUV_factors, model_str_full,model_str_red, opt_DESeq2_parallel,min_p_donors_groups){

  ## generate indices for random dataset halves
  idx_rnd <- get_indices_for_random_dataset_halves(DDS,min_p_donors_groups)
  
  #initialize correlation value
  M_cor <- matrix(NaN,1,2)
  
  if (n_RUV_factors>0){
    ## (1) calculate RUVs for n_RUV_factors factors based on filtered data set (including all samples) 
    #set <- RUVg(D$set_ini, D$empirical, k=n_RUV_factors)
    
    #load RUVs
    set <- readRDS(file=paste0(paths$results_RUV,"set_",as.character(n_RUV_factors),"_RUVS.rds"))
    
    RUV_factors <- rep(NaN,n_RUV_factors)
    for (k in seq(1,n_RUV_factors)){
      RUV_factors[k] <- paste0("W_",as.character(k))
    }
    
    ## (2) update full and reduced design, include n_RUV_factors RUV factors:
    model_str_full <- get_design_str(factors_to_include=RUV_factors,opt_RUV=TRUE,opt_red=FALSE)
    model_str_red <- get_design_str(factors_to_include=RUV_factors,opt_RUV=TRUE,opt_red=TRUE)
  }
  
  design_formula_full <- formula(eval(parse(text=model_str_full)))
  design_formula_red <- formula(eval(parse(text=model_str_red)))
  
  #initialize data halve DDS objects
  DDS_1st_half_k <- c()
  DDS_2nd_half_k <- c()
  
  #split DDS object and update design:
  DDS_1st_half_k <- DDS[,idx_rnd$first]
  DDS_2nd_half_k <- DDS[,idx_rnd$second]
  DDS_1st_half_k@design <- design_formula_full
  DDS_2nd_half_k@design <- design_formula_full
  
  #for each half: determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
  sfe_type_1st <- get_sfe_type(DDS_1st_half_k)
  sfe_type_2nd <- get_sfe_type(DDS_2nd_half_k)
  
  print("Check 2!")
  
  if (n_RUV_factors != 0){
    # integrate RUV factors in DDS objects colData
    DDS_1st_half_k <- integrate_RUV_factors_in_DDS_obj(DDS_1st_half_k,set,idx_rnd,'first')
    DDS_2nd_half_k <- integrate_RUV_factors_in_DDS_obj(DDS_2nd_half_k,set,idx_rnd,'second') 
  }
  
  print("Check 3!")
  
  flag_1st <- tryCatch(
    {
      #run DSEq2:
      DDS_1st_half_k <- DESeq(object=DDS_1st_half_k,
                              parallel=opt_DESeq2_parallel,
                              sfType = sfe_type_1st, 
                              test="LRT", 
                              reduced = design_formula_red, 
                              useT=TRUE, 
                              minmu = 1e-6, 
                              minReplicatesForReplace = Inf)
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
      flag <- 0
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      flag <- 303
      return(flag)
    },
    warning=function(cond) {
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      flag <- 101
      return(flag)
    },
    finally={
      #store correlation of log2FoldChanges in correlation matrix (n_splits x n_RUV_factors_max)
      message(paste0("\n n_RUV_factors = ",as.character(n_RUV_factors))," is done (1/2)! \n")
    }
  )
  print("Check 4!")
  flag_2nd <- tryCatch(
    {
      #run DSEq2:
      DDS_2nd_half_k <- DESeq(object=DDS_2nd_half_k,
                              parallel=opt_DESeq2_parallel,
                              sfType = sfe_type_2nd, 
                              test="LRT", 
                              reduced = design_formula_red, 
                              useT=TRUE, 
                              minmu = 1e-6, 
                              minReplicatesForReplace = Inf)
      
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped inside a function (unlike that
      # for the condition handlers for warnings and error below)
      flag <- 0
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      flag <- 303
      return(flag)
    },
    warning=function(cond) {
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      flag <- 101
      return(flag)
    },
    finally={
      #store correlation of log2FoldChanges in correlation matrix (n_splits x n_RUV_factors_max)
      message(paste0("\n n_RUV_factors = ",as.character(n_RUV_factors))," is done (2/2)! \n")
    }
  )
  print("Check 5!")
  alpha_val<-0.05 # doesn't really matter since all log2FC are used for correlation calculation anyways
  if (flag_1st==0 & flag_2nd==0){
    res_1st_half_k <- results(DDS_1st_half_k,alpha=alpha_val,contrast=c("Group", "SCZ","CTRL"),independentFiltering = FALSE)
    res_2nd_half_k <- results(DDS_2nd_half_k,alpha=alpha_val,contrast=c("Group", "SCZ", "CTRL"),independentFiltering = FALSE)
    # calculate spearman/ pearson correlation
    M_cor[1,1] <- cor(res_1st_half_k$log2FoldChange,res_2nd_half_k$log2FoldChange,method="spearman")
    M_cor[1,2] <- cor(res_1st_half_k$log2FoldChange,res_2nd_half_k$log2FoldChange,method="pearson")
  }
  
  return(M_cor)
}

#get background gene set
get_background <- function(path_background,BrainCortex_TH,genetype,opt_background_data,opt_filter_for_protein_coding_genes,version_nr){
  if (version_nr=='v3'){
    gene_matrix <- LOAD("geneMatrix_v2.tsv",path_background)
  }else{
    gene_matrix <- LOAD(paste0(paste0("geneMatrix_",version_nr),".tsv"),path_background)
  }
  if (opt_background_data=="cortex"){
    if (opt_filter_for_protein_coding_genes==TRUE){
      b <- gene_matrix %>% filter(gene_type == genetype, #only protein coding genes
                                  BrainCortex>=BrainCortex_TH,
                                  !is.na(BrainCortex),
                                  PAR == FALSE,
                                  dup_gene_name == FALSE, # only non-duplicate gene names
                                  hg38h0 != "chrM")#, # remove if mapped to mitochondrial chromosome??
    } else {
      b <- gene_matrix %>% filter(BrainCortex>=BrainCortex_TH,
                                  !is.na(BrainCortex),
                                  PAR == FALSE,
                                  dup_gene_name == FALSE, # only non-duplicate gene names
                                  hg38h0 != "chrM")#, # remove if mapped to mitochondrial chromosome??
    }
  } else if (opt_background_data=="SnRNAseqDet"){
    genes_snRNAseq_detected <- LOAD("background_genes_from_snRNAseq.csv",path_background)
    if (opt_filter_for_protein_coding_genes==TRUE){
      protein_coding <- gene_matrix %>% filter(gene_type == genetype, #only protein coding genes
                                               PAR == FALSE,
                                               dup_gene_name == FALSE, # only non-duplicate gene names
                                               hg38h0 != "chrM")#, # remove if mapped to mitochondrial chromosome??
      b <- genes_snRNAseq_detected %>% filter(ensgid %in% protein_coding$ensgid)
    } else{
      b <- genes_snRNAseq_detected
    }
  } else if (opt_background_data=="genesMappingToProteinsDetectedInProteomics"){
    b <- LOAD("background_genes_from_proteomics.csv",path_background)
  }
  return(b)
}

calculate_and_export_RUVs = function(D,DDS_RUV,colData,n_RUV_factors,path_res_RUV,opt_plot_variance, opt_downsampled_data,i_ds){

  set_ini <- RUVg(D$set_ini, D$empirical, k=n_RUV_factors)
  set <- betweenLaneNormalization(set_ini, which="upper")
  
  #obtain and export normalized counts:
  RUV_normalized_counts <- normCounts(set)
  setwd(path_res_RUV)
  
  if (opt_downsampled_data==TRUE){
    file_ending <- paste0("_",as.character(i_ds),".Rds")
  } else{
    file_ending <- ".Rds"
  }
  saveRDS(RUV_normalized_counts,file=paste0("Counts_RUV_normalized_",n_RUV_factors,"_RUV_factors",file_ending))
  
  #perform RUV analysis to test if design in deseq2 should be more complicated or if simple design is sufficient
  #Removal of unwanted variation
  #pull out a set of empirical control genes 
  #genes that do not have a small p-value in initial DESeq2 run and low CT specificity
  
  #plot RLE = log-ratio of read count to median read count across sample
  plot_RUV_RLE(set,DDS_RUV,path_res_RUV,n_RUV_factors)
  
  plot_RUV_PCA(set,DDS_RUV,path_res_RUV,n_RUV_factors)
  
  #there should be no differences between the groups
  # test if differences in factor between groups (t-test) and exclude any factor with corrected p-val <0.05
  factors_to_remove <- c()
  print("Raw p-values for t-test with H0: RUV factor i does not differ in mean between the two groups SCZ and CTRL: \n")
  for (i in seq(1,n_RUV_factors)){
    #plot factors estimated by RUV:
    plot_RUV_factors(DDS_RUV,pData(set),data_all$celltype,path_res_RUV,'Group',n_RUV_factors,i,opt_downsampled_data,i_ds)
    
    factor_str <- paste0("W",i)
    eval(parse(text=paste0(factor_str,"<-as.numeric(pData(set)[,",as.character(i),"])")))
    test_res <- t.test(eval(parse(text=factor_str))[DDS_RUV$Group=='CTRL'],eval(parse(text=factor_str))[ DDS_RUV$Group=='SCZ'],alternative = "two.sided", 
                       mu = 0, paired = FALSE, var.equal=FALSE, conf.level = 0.95)
    print(test_res$p.value)
    if (test_res$p.value<0.05/n_RUV_factors){
      factors_to_remove <- c(factors_to_remove,factor_str)
    }
  }
  #analyse each factors contribution to the expression variation of each gene 
  #store RUV factors in dataframe

  a <- seq(1,n_RUV_factors)
  RUV_factors <- eval(parse(text=paste0('data.frame(',paste(paste0("W",a),collapse=','),')')))
  rownames(RUV_factors)<-colData$donor_ID_python
  #rownames(RUV_factor_info) <- 
  #expression values can be accessed via: counts(DDS)
  #specify Wi as fixed effects:
  #form <- ~ W1 + W2 + W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10
  form <- paste0('~',paste(paste0("W",a),collapse='+'))
  
  if (opt_plot_variance==TRUE){ ##much slower as running without
    #continue here: which data to provide
    #variance stabilizing transformation
    nsub = sum( rowMeans( counts(DDS_RUV, normalized=TRUE)) > 5 )
    
    counts_vst <- vst(counts(set),nsub = nsub)
    
    varPart <- fitExtractVarPartModel(counts_vst, form, RUV_factors)
    
    # sort variables (i.e. columns) by median fraction
    # of variance explained
    vp <- sortCols( varPart )
    
    # Bar plot of variance fractions for the first 10 genes
    plot_RUV_var_fractions_per_gene(vp,path_res_RUV,n_RUV_factors,opt_downsampled_data,i_ds)
    
    # violin plot of contribution of each variable to total variance
    plot_RUV_var_per_factor(vp,path_res_RUV,n_RUV_factors,opt_downsampled_data,i_ds)
    
    #cumulative variance:
    plot_RUV_cumulative_var(vp,path_res_RUV,n_RUV_factors,opt_downsampled_data,i_ds)
  }
  
  #if one wants to use them: store result and load for DEG analysis
  #save RUV factors
  setwd(path_res_RUV)
  saveRDS(RUV_factors,file=paste0("RUV_factors_",n_RUV_factors,"_factors",file_ending))
  
  return(RUV_factors)
}

get_design_str = function(factors_to_include,opt_RUV,opt_red){
  design_str <- "~"
  if (opt_RUV==TRUE){
    #if (opt_red==FALSE){
    #  design_str <- paste0("design(",DDS_obj_str,") <- formula( ",design_str)
    #}
    for (wi in factors_to_include){
      #add the wi as variables to DDS object: 
      #DDS_RUV$W1 <- RUV_factor_info$W1
      #eval(parse(text=paste0(DDS_obj_str,"$",wi," <- RUV_factor_info$",wi)))
      #build design formula
      #design(DDS_RUV) <- formula(~ W1 + W2 + ... + Group)
      if (opt_red==TRUE){
        if (wi==factors_to_include[length(factors_to_include)]){
          design_str<- paste0(design_str,wi)
        } else{
          design_str <- paste0(design_str,wi," + ")
        }
      } else{
        design_str <- paste0(design_str,wi," + ")
      }
    }
  } else{
    if (opt_red==TRUE){
      design_str <- design_str
    } else{
      design_str <- paste0("design(",DDS_obj_str,") <- formula( ",design_str)
    }
  }
  
  if (opt_red==FALSE){
    design_str <- paste0(design_str," Group")
  } else{
    if (design_str == "~"){
      design_str <- "~ 1"
    }
  }
  
  return(design_str)
}

# integrate RUV factors in DDS objects colData
integrate_RUV_factors_in_DDS_obj = function(DDS_RUV,set,idx_rnd,dataset_half){
  add_RUV_str_colData <- ''
  DDS_RUV_names <- names(DDS_RUV@colData@listData)
  if  (dataset_half=='none'){
    RUV_names <- colnames(set)
    str1 = ", "
    str2 = ""
  }else if (dataset_half=='first'){
    RUV_names <- colnames(pData(set))
    str1="[idx_rnd$first], "
    str2="[idx_rnd$first]"
  }else if (dataset_half=='second'){
    RUV_names <- colnames(pData(set))
    str1 = "[idx_rnd$second], "
    str2 = "[idx_rnd$second]"
  }
  
  for (RUV_i in RUV_names){
    if (RUV_i != RUV_names[length(RUV_names)]){
      add_RUV_str_colData <- paste0(add_RUV_str_colData,RUV_i,'=set$',RUV_i,str1)
    } else{
      add_RUV_str_colData <- paste0(add_RUV_str_colData,RUV_i,'=set$',RUV_i,str2)
    }
    DDS_RUV_names <- c(DDS_RUV_names,RUV_i)
  }
  DDS_RUV@colData@listData <- c(DDS_RUV@colData@listData,
                                eval(parse(text=paste0(paste0('list(',add_RUV_str_colData),')'))))
  names(DDS_RUV@colData@listData) <- DDS_RUV_names
  
  return(DDS_RUV)
}

get_settings_for_input_data<- function(opt_input_data){
  if (opt_input_data == "Patient_clustering_genes_GSA_result"){
    n_clusters_range = c('SCZ_cases')#,'whole_dataset','controls')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('cortex/','SnRNAseqDet/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c('all')
  } else if (opt_input_data =="Patient_clustering_genes_GSEA_result"){
    n_clusters_range = c('controls','SCZ_cases','whole_dataset')
    min_overlap_range <- c(0)
    background_folders <- c('none/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c('all')
  } else if (opt_input_data == "DEG_GSA_result"){ 
    n_clusters_range <- c(15)#,37,3)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('cortex/','cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c("Excitatory_Layer_5_6_IT_neurons_I","Excitatory_Layer_5_6_IT_neurons_II",
                   "Excitatory_Layer_5_6_CT_and_NP_neurons","Excitatory_Layer_2_3_IT_neurons_I",
                   "Excitatory_Layer_2_3_IT_neurons_II", "Excitatory_Layer_3_4_IT_neurons","Excitatory_Layer_3_6_IT_neurons",
                   "Astrocytes","Endothelial_and_mural_cells","Microglial_cells","Oligodendrocyte_progenitor_cells","Oligodendrocytes",
                   "Inhibitory_LAMP5_neurons","Inhibitory_PVALB_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons")                     
    #old definition:                                        
    #CT_groups <- c('all','Astrocytes', "Microglial_cells","Oligodendrocyte_progenitor_cells",'Neurons')
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){ 
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('norm_basic/')
    min_overlap_range <- c(0)
    background_folders <- c('none/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- TRUE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
    n_clusters_range <- c(15)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('cortex/')#,'cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if(opt_input_data == "all_GRN_modules_gprofiler_result"|| opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    n_clusters_range <- c(15)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('cortex/')#,'cortex_pc/','SnRNAseqDet/','SnRNAseqDet_pc/')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "all_GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('GSA/')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- TRUE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if (opt_input_data == "GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('GSA/')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- FALSE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data =="longread_results"){
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "TWAS_result"){
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "proteomics_result"){
    background_folders <- c('cortex/')
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "DEG_GSA_result_all_CTs"){
    n_clusters_range <- c(15)#,37,3)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('cortex/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c("Excitatory_Layer_5_6_IT_neurons_I","Excitatory_Layer_5_6_IT_neurons_II",
                   "Excitatory_Layer_5_6_CT_and_NP_neurons","Excitatory_Layer_2_3_IT_neurons_I",
                   "Excitatory_Layer_2_3_IT_neurons_II", "Excitatory_Layer_3_4_IT_neurons","Excitatory_Layer_3_6_IT_neurons",
                   "Astrocytes","Endothelial_and_mural_cells","Microglial_cells","Oligodendrocyte_progenitor_cells","Oligodendrocytes",
                   "Inhibitory_LAMP5_neurons","Inhibitory_PVALB_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons")
  }
  if ((str_detect(opt_input_data,"GRN")) && (opt_input_data != "all_GRN_modules_gprofiler_result") && (opt_input_data != "all_bordeaux_modules_gprofiler_result")){
    n_clusters_range <- get_cell_type_clusters_GRN(opt_metacell_settings)
  }
  return(list(n_clusters_range,sub_folders,min_overlap_range,background_folders,opt_sim_based_on,CT_groups))
}

get_GSA_results_as_DF <- function(data_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data){
  setwd(data_path)
  D <- c()
  if (opt_up_and_down_sep==FALSE){
    bool_first = TRUE
    file_list <- Sys.glob("*__all.csv")
    # opt_input_data=="GRN_modules_GSA_result" | opt_input_data=="all_GRN_modules_GSA_result")
    for (file in file_list){
      d <- read.csv2(file) 
      d <- d %>% select(c('group', 'subgroup', 'geneset','P.bonf.group','TestVar')) %>% filter(P.bonf.group<=p_val_cutoff)
      if (bool_first==FALSE){
        D <- rbind(D,d)
      } else{
        D <- d
        bool_first <- FALSE
      }
    }
  } else{
    D$up <- c()
    D$down <- c()
    bool_first = TRUE
    reg_strings <- c('up','down')
    for (reg_str in reg_strings){
      file_list <- Sys.glob(paste0("*",reg_str,"_all.csv"))
      for (file in file_list){
        d <- read.csv2(file) 
        d <- d %>% select(c('group', 'subgroup', 'geneset','P.bonf.group','TestVar')) %>% filter(P.bonf.group<=p_val_cutoff)
        if (bool_first==FALSE){
          if (reg_str == 'up'){
            D$up <- rbind(D$up,d)
          } else{
            D$down <- rbind(D$down,d)
          }
        } else{
          if (reg_str == 'up'){
            D$up <- d
          } else{
            D$down <- d
          }
          bool_first <- FALSE
        }
      }
    }
  }
  return(D)
}

get_gprofiler_results_as_DF <- function(data_path, p_val_cutoff){
  file <- Sys.glob(paste0(data_path,"g_profiler_results_*.csv"))
  d <- read.csv(file) 
  module_columns <- grep('_module_',colnames(d),value=TRUE)
  module_columns_short <- module_columns[endsWith(module_columns,'_grey')==FALSE]
  module_columns_super_short <- as.character(lapply(strsplit(module_columns_short, split="_"),tail, n=1))
  d <- d %>% select(c('source','term_id', 'term_name', module_columns_short))
  colnames(d) <- c('source','term_id', 'term_name', module_columns_super_short)
  #filter modules for SCZ sign? 
  scz_sign<- read.csv(paste0(data_path,"DMEs_LR.csv"))
  sign_modules <- scz_sign$module[scz_sign$p_val_adj<=0.05]
  d <- d %>% select(c('source','term_id', 'term_name', sign_modules))
  D <- pivot_longer(d, cols= sign_modules)
  #rename columns
  colnames(D)[colnames(D)=="name"] <- "Module"
  colnames(D)[colnames(D)=="value"] <- "pvalue"
  colnames(D)[colnames(D)=="term_id"] <- "ID"
  colnames(D)[colnames(D)=="source"] <- "subgroup"
  #filter for 
  D <- D %>% filter(pvalue<=p_val_cutoff)
  return(D)
}

get_longread_results <- function(data_path, p_val_cutoff){
  filename <- paste0(data_path,"go_hypergeoPathway_filtE_fdr10.tsv")
  d <- read.table(file=filename, sep = '\t', header = TRUE)
  d <- d %>% select(c('subgroup', 'geneset','Phyper',"BH")) 
  #d_GO_BP <- dplyr::filter(d, grepl("GO:BP",subgroup))
  #d_GO_CC <- dplyr::filter(d, grepl("GO:CC",subgroup))
  #d_GO_MF <- dplyr::filter(d, grepl("GO:MF",subgroup))
  #d_GO_BP_sign <- d_GO_BP %>% filter(Phyper/dim(d_GO_BP)[1] <= p_val_cutoff)
  d$P.bonf.group <- d$BH
  d$ID <- sub(" .*", "", d$geneset)
  d$TestVar <- "longread"
  d <- d[ , -which(names(d) %in% c("Phyper","BH","group"))]
  return(d)
}

get_TWAS_results <- function(data_and_results_path, p_val_cutoff){
  filename <- paste0(data_and_results_path, "df_15CT_TWAS0.3_GSA0.05_GO_group3class.tsv")
  d <- read.table(file=filename, sep = '\t', header = TRUE)
  d <- d %>% select(c('subgroup', 'geneset','P.fdr.group')) 
  d$P.bonf.group <- d$P.fdr.group
  d$ID <- sub(" .*", "", d$geneset)
  d$TestVar <- "TWAS"
  d <- d[ , -which(names(d) %in% c("P.fdr.group"))]
  return(d)
}

get_proteomics_results <- function(data_and_results_path, p_val_cutoff){
  filename <- paste0(data_and_results_path, "proteomics_both_sign.csv")
  d <- read.table(file=filename, sep = ';', header = TRUE)
  d$ID <- sub(" .*", "", d$geneset)
  d$P.bonf.group <- d$P.fdr.group
  d <- d %>% select(c('subgroup', 'geneset','ID','P.bonf.group')) 
  d$TestVar <- "proteomics"
  return(d)
}

get_data <- function(opt_input_data,main_path, data_and_results_path,p_val_cutoff, cf, sf, opt_up_and_down_sep,opt_metacell_settings){
  #"DEG_GSA_result_all_CTs"
  ##"longread_results"
  #"TWAS_result"
  ##"all_GRN_modules_gprofiler_result"
  #"all_bordeaux_modules_gprofiler_result"
  #"proteomics_result"
  
  if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    
    if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
      D_d <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
    }
    sf <- 'norm_basic/'
    n_c_r <- get_cell_type_clusters_GRN(opt_metacell_settings)
    
    bool_first <- TRUE
    for (n_c in n_c_r){
      cf <- paste0(n_c,'/')
      data_and_results_path_GRN <- paste0(data_and_results_path,cf,sf)
      D_m_tibble_tmp <- get_gprofiler_results_as_DF(data_and_results_path_GRN, p_val_cutoff)
      D_m_tibble_tmp <- D_m_tibble_tmp %>% mutate_at(.vars = "Module", ~ paste0(n_c,"_",.x))
      D_m_tibble_tmp <- D_m_tibble_tmp %>% unite("geneset", ID:term_name, sep=" ", na.rm=TRUE)
      if (bool_first==TRUE){
        D_m_tibble <- D_m_tibble_tmp
        bool_first <- FALSE
      } else{
        D_m_tibble <- rbind(D_m_tibble,D_m_tibble_tmp)
      }
    }
    #harmonize dataframes:
    D_m <- as.data.frame(D_m_tibble)
    if (nrow(D_m)){
      colnames(D_m)[colnames(D_m)=="Module"] <- "TestVar"
      D <- D_m
    }
    if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
      if (nrow(D_d$up)>0){
        colnames(D_d$up)[colnames(D_d$up)=="P.bonf.group"] <- "pvalue"
        D_d$up <- subset(D_d$up,select= -c(group))
        if (nrow(D)>0){
          D <- rbind(D,D_d$up)
        } else{
          D <- D_d$up
        }
      }
      if (nrow(D_d$down)>0){
        colnames(D_d$down)[colnames(D_d$down)=="P.bonf.group"] <- "pvalue"
        D_d$down <- subset(D_d$down,select= -c(group))
        if (nrow(D)>0){
          D <- rbind(D,D_d$down)
        } else{
          D <- D_d$down
        }
      }
    } 
    if (opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "all_bordeaux_modules_gprofiler_result"){
      
      if (opt_input_data=="all_bordeaux_modules_gprofiler_result"){
        #filter for bordeaux modules:
        D <- D[which(D$TestVar %in% c("Oligodendrocytes_red","Excitatory_Layer_5-6_IT_neurons_I_green",
                                      "Inhibitory_LAMP5_neurons_red","Excitatory_Layer_2-3_IT_neurons_II_black", 
                                      "Astrocytes_yellow","Inhibitory_PVALB_neurons_magenta",
                                      "Excitatory_Layer_2-3_IT_neurons_I_red","Inhibitory_VIP_neurons_pink")), ]
        D$TestVar <- "bordeaux_modules"
      } else{
        D$TestVar <- "modules"
      }
      #D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
      D$P.bonf.group <- D$pvalue
      D <- D[ ,-which(names(D) %in% c("pvalue"))]
      D$ID <- sub(" .*", "", D$geneset)
    }
    
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){ # implementation for gprofiler output
    D <- get_gprofiler_results_as_DF(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "longread_results"){
    D <- get_longread_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "TWAS_result"){
    D <- get_TWAS_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "proteomics_result"){
    D <- get_proteomics_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "GRN_modules_GSA_result"){
    D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
  } else{ # implementation for GSA output
    D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
    if (opt_input_data=="DEG_GSA_result_all_CTs"){
      if (opt_up_and_down_sep==TRUE){
        D$up$TestVar <- "DEGs_up"
        D$down$TestVar <- "DEGs_down"
        D <- rbind(D$up,D$down)
      } else{
        D$TestVar <- "DEGs"
      }
      D <- D[ ,-which(names(D) %in% c("group"))]
      #add GO ID column:
      D$ID <- sub(" .*", "", D$geneset)
    }
  }
  return(D)
}

get_all_data <- function(opt_input_data_vec,opt_metacell_settings, main_path,p_val_cutoff){
  bool_first <- TRUE
  for (opt_input_data in opt_input_data_vec){
    #print(opt_input_data)
    opt <- get_settings_for_input_data(opt_input_data)
    n_clusters_range <- opt[[1]]
    sub_folders <- opt[[2]]
    min_overlap_range <- opt[[3]]
    background_folders <- opt[[4]]
    opt_sim_based_on <- opt[[5]]
    CT_groups <- opt[[6]]
    #"DEG_GSA_result_all_CTs"
    #"longread_results"
    #"TWAS_result"
    #"all_GRN_modules_gprofiler_result"
    #"all_bordeaux_modules_gprofiler_result"
    #"proteomics_result"
    if (opt_input_data == "GRN_modules_gprofiler_result" | opt_input_data =="all_GRN_modules_gprofiler_result" | opt_input_data =="all_bordeaux_modules_gprofiler_result" | opt_input_data == "all_GRN_modules_GSA_result" | opt_input_data == "GRN_modules_GSA_result"){ 
      if (opt_with_MT_genes){
        MT_subfolder  <- c(paste0(opt_metacell_settings,"/with_MT_genes/"))
      } else{
        MT_subfolder <- c(paste0(opt_metacell_settings,"/without_MT_genes/"))
        if (opt_new_profiler_settings==TRUE){
          MT_subfolder <- c(paste0(opt_metacell_settings,"/without_MT_genes_new_gprofiler_settings/"))
        }
      }
    } else{
      MT_subfolder <- ""
    }
    if (opt_input_data=="GRN_modules_gprofiler_result"){
      opt_up_and_down_sep = FALSE #
    } else if (opt_input_data=="Patient_clustering_genes_GSA_result" | opt_input_data=="GRN_modules_GSA_result" | opt_input_data=="all_GRN_modules_GSA_result" | opt_input_data =="all_bordeaux_modules_gprofiler_result" ){
      opt_up_and_down_sep = FALSE #
    } else{
      opt_up_and_down_sep = TRUE#FALSE #
    }
    for (n_clusters in n_clusters_range){
      if (opt_input_data %in% c("Patient_clustering_genes_GSA_result","Patient_clustering_genes_GSEA_result")){
        cluster_folder <- paste0(n_clusters,'/no_normalization/batch_correction_for_Library/')#norm_var_p_vals_01_var_scale_power_1/') #norm_var_p_vals_00001/')#
        #if (n_clusters=='controls'){
        sub_folders <- c('factor_1/','factor_2/','factor_3/')
        #} else{
        #  sub_folders <- c('factor_1/','factor_2/','factor_3/','factor_4/','factor_5/','factor_6/')
        #}
      } else if (opt_input_data == "DEG_GSA_result" || opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data =="all_bordeaux_modules_gprofiler_result" || opt_input_data == "DEG_GSA_result_all_CTs"){  
        if (n_clusters==3){
          cluster_folder <- paste0(n_clusters,'_classes/')
        } else{
          cluster_folder <- paste0(n_clusters,'_CTs/')
        }
      } else if (opt_input_data=="GRN_modules_gprofiler_result"){
        cluster_folder <- paste0(n_clusters,'/')
      } else{
        cluster_folder <- ""
      }
      for (sub_folder in sub_folders){
        for (min_overlap in min_overlap_range){
          for (background_folder in background_folders){
            
            #load significant pathways:
            data_and_results_path <- get_results_path(main_path,opt_input_data,cluster_folder,min_overlap,background_folder,sub_folder,opt_metacell_settings)
            D <- get_data(opt_input_data,main_path, data_and_results_path,p_val_cutoff, cluster_folder, sub_folder, opt_up_and_down_sep,opt_metacell_settings)
            #make sure order of columns is harmonized:
            D <- D[ ,c("subgroup","geneset","ID","P.bonf.group","TestVar")]
            if (bool_first==TRUE){
              D_all <- D
              bool_first <- FALSE
            } else{
              D_all <- rbind(D_all,D)
            }
          }
        }
      }
    }
  }
  return(D_all)
}
get_results_path <- function(main_path,opt_input_data,cluster_folder,min_overlap,background_folder,sub_folder,opt_metacell_settings){
  if (opt_input_data == "DEG_GSA_result"|| opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "DEG_GSA_result_all_CTs"){
    data_and_results_path <- paste0(main_path,'/output/GSA_analysis/DEGs/v3/',cluster_folder,sub_folder,'min_overlap_',as.character(min_overlap),'/6_RUVs/',background_folder)
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){
    data_and_results_path <- paste0(main_path,'/output/cluster_name_15CTs/scz/cor_bicor/',cluster_folder,sub_folder)
  } else if (opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    data_and_results_path <- paste0(main_path,'/output/cluster_name_15CTs/scz/cor_bicor/',cluster_folder,"/")
  } else if(opt_input_data =="longread_results"){
    data_and_results_path <- paste0(main_path,'/output/GSA_analysis/DIUGs/v3/alpha3/min_overlap_3/',background_folder,"/")
  } else if (opt_input_data == "TWAS_result"){
    data_and_results_path <- paste0(main_path,"/output/TWAS_pathways/")
  } else if (opt_input_data == "proteomics_result"){
    data_and_results_path <- paste0(main_path,"/output/GSA_analysis/DAPs/v3/alpha3/min_overlap_3/",background_folder,"/")
  } else if (opt_input_data== "GRN_modules_GSA_result"){
    data_and_results_path <- paste0(main_path,'/output/cluster_name_15CTs/scz/cor_bicor/GSA_result/',cluster_folder,"/",sub_folder)
  } else if (opt_input_data == "all_GRN_modules_GSA_result"){
    data_and_results_path <- paste0(main_path,'/output/cluster_name_15CTs/scz/cor_bicor/GSA_result/',sub_folder)
  }
  #print(data_and_results_path)
  return(data_and_results_path)
}


plot_pathway_clustering <- function(simMatrix, reducedTerms, results_path, add_str, GO_category, reg_str, ct_group,opt_up_and_down_sep,opt_input_data){
  if (opt_input_data != "GRN_modules_gprofiler_result"){
    if (opt_up_and_down_sep==TRUE){
      subfolder = 'up_and_down_sep/'
    }else{
      subfolder = 'up_and_down_together/'
    }
  } else{
    subfolder = ''
  }
  path <- paste0(results_path,"figures/rrvgo/",subfolder)
  #make sure path exists:
  dir.create(path, recursive=TRUE, showWarnings = FALSE)
  
  graphic_filename_1 <- paste0(path,"heatmap_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_1,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(heatmapPlot(simMatrix,
                    reducedTerms,
                    annotateParent=TRUE,
                    annotationLabel="parentTerm",
                    fontsize=6))
  dev.off()
  
  if (dim(reducedTerms)[1]>4){
    graphic_filename_2 <- paste0(path,"scatter_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
    pdf(file = graphic_filename_2,   # The directory you want to save the file in
        width = 6, # The width of the plot in inches
        height = 6) # The height of the plot in inches
    print(scatterPlot(simMatrix, reducedTerms, addLabel=TRUE))
    dev.off()
  }
  
  
  # module_list <- treemapPlot(reducedTerms,size="score")
  # #change some colors
  # if (ct_group=="all" & reg_str=='down' & GO_category=='BP'){
  #   module_list <- hand_curate_colors(module_list)
  #   browser()
  # }
  graphic_filename_3 <- paste0(path,"treemap_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_3,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(module_list <- treemapPlot(reducedTerms,size="score"))
  dev.off()
  
  graphic_filename_4 <- paste0(path,"wordcloud_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_4,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(wordcloudPlot(reducedTerms, min.freq=1, colors="black"))
  dev.off()
  
  return(module_list)
}


get_cell_type_clusters_GRN <- function(opt_metacell_settings){
  if (opt_metacell_settings == "50_13_250" ){
    n_c_r = c("Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
              "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
              "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
              "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I")  
  } else if (opt_metacell_settings == "30_6_120"){
    n_c_r = c("Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
              "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
              "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
              "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I",
              "Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells")  
  }
  return(n_c_r)
}

#plot and save significant genes of different settings 
plot_significant_genes <- function(DDS,res,alpha,mode,path_results,i) {
  sign_genes=which(res$padj<alpha)
  if (sum(res$padj < alpha_val, na.rm=TRUE)>0){
    #n = ceiling(sqrt(length(sign_genes)))
    filename <- paste0(path_results,"/Sign_genes_",mode,sep='_',".png")
    png(file=filename,width=600, height=550, units='mm', res=300)
    par(mfrow=c(1,1),mar=c(1,1,1,1)) 
    for (sg in sign_genes){
      plotCounts(DDS, gene=sg, intgroup="Group",cex=2, cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    }
    dev.off()
    print("significant genes:")
    print(rownames(res)[sign_genes])
  }
}

#plot size factors sorted by group
plot_size_factors <- function(DDS,celltype,path_results,sfe_type,i) {
  filename <- paste0(path_results,"Size_factors.png")
  png(file=filename,width=600, height=550, units='mm', res=300)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  boxplot(DDS$sizeFactor ~ DDS$Group, notch=F, col = c("steelblue","orange"),main=paste0(paste0(paste0(sfe_type," ")," size factors for "),celltype))
  stripchart(DDS$sizeFactor[DDS$Group=="CTRL"] ~ DDS$Group[DDS$Group=="CTRL"], add=T, pch=21, vertical = TRUE, method='jitter', jitter=0.2, col = c("black"),bg = c("steelblue"))
  stripchart(DDS$sizeFactor[DDS$Group=="SCZ"] ~ DDS$Group[DDS$Group=="SCZ"], add=T, pch=21, vertical = TRUE, method='jitter', jitter=0.2, col = c("black"),bg = c("orange"))
  dev.off()
  graphics.off()
}

#plot RUV analysis factors sorted by group
plot_RUV_factors <- function(DDS,pData_set,celltype,path_results,variable_name,n_RUV_factors,i,opt_downsampled_data,i_ds) {
  if (opt_downsampled_data==TRUE){
    file_ending <- paste0(i_ds,".pdf")
  } else{
    file_ending <- ".pdf"
  }
  filename <- paste0(path_results,"/RUV_factors_",variable_name,n_RUV_factors,file_ending)
  pdf(file=filename)
  if (dim(pData_set)[2]>4){
    par(mfrow=c(2,as.integer(ceiling(dim(pData_set)[2])/2)),mar=c(1,1,1,1)+3.5)
  }
  else{
    par(mfrow=c(1,dim(pData_set)[2]),mar=c(1,1,1,1)+3.5)
  }
  cont_vars <- c('Age','PMI','mean_reads_per_umi','num_umi','p_cell_del_filt','num_reads')
  cat_vars <-c('Group','Sex','Library')
  #par(mfrow = c(2, 1), mar = c(3,5,3,1))
  
  #categorical variables:
  for (i in seq(1,dim(pData_set)[2])){
    if (variable_name %in% cat_vars){
      if (variable_name=='Group'){
        colors <- c("steelblue","orange")
        #boxplot(pData_set[, i] ~ DDS$Group, notch=F, col = colors, main=paste0(paste0(paste0("W",i)," for "),celltype),las=2,outline=FALSE,xlab="Group",ylab=paste0("W",i))
        #stripchart(pData_set[,i] ~ DDS$Group, group.names = unique(DDS$Group),add=T,pch=21, vertical = TRUE,method='jitter', jitter=0.2, col="black", bg = "white", xlab="Group",ylab=paste0("W",i))
      }
      if (variable_name=='Sex'){
        colors = c("cornflowerblue","deeppink3")
      }
      if (variable_name == 'Library'){
        colors <- c("chocolate1","firebrick3","aquamarine3","dodgerblue4","darkorchid")
      }
      x <- eval(parse(text=paste0("DDS$",variable_name)))
      boxplot(pData_set[, i] ~ x, notch=F, col = colors, main=paste0(paste0("W",i)," for "),las=2,outline=FALSE,xlab=variable_name,ylab=paste0("W",i))
      stripchart(pData_set[,i] ~ x, group.names = unique(x),add=T,pch=21, vertical = TRUE,method='jitter', jitter=0.2, col="black", bg = "white", xlab=variable_name,ylab=paste0("W",i))
    }
    
    #continuous variables:
    if (variable_name %in% cont_vars){
      if (variable_name=='p_cell_del_filt'){
        x <- DDS$p_cell_del_filt[1:82]
      } else{
        x <- eval(parse(text=paste0("DDS$",variable_name)))
      }
      y <- pData_set[, i]
      rho_x_y <- round(cor(x,y,use="pairwise.complete.obs",method='spearman'),2)
      # Change point shape (pch = 19) and remove frame.
      plot(x, y, main = paste0("rho=",rho_x_y),xlab = variable_name, ylab = paste0("W",i))
    }
  }
  dev.off()
  graphics.off()
}

#MA-plot: log2 fold changes attributable to a given variable (Group) over the mean of normalized counts for all the samples in the DESeqDataSet
#         It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes 
#         from low count genes without requiring arbitrary filtering thresholds.
plot_MA_for_various_shrinkage_types <- function(RES_obj,shrinkage_types,path_res){
  png(file=paste0(path_res,"Mean_norm_counts_vs_LFC.png"),width=600, height=200, units='mm', res=300)
  par(mfrow=c(1,length(shrinkage_types)),mar=c(1,1,1,1)) 
  for (st in shrinkage_types){
    if (st == "no_shrinkage"){
      main_str="No shrinkage of effect size"
    } else{
      main_str=paste0("Shrinkage with ",st)
    }
    plotMA(RES_obj[[st]], ylim=c(-2,2), main=main_str)#
  }
  dev.off()
}

#plot boxplots and correlation plot of metadata:
visualize_metadata <- function(metadata,path_results){
  #for numeric columns plot correlation:
  filename <- paste0(paste0(path_results,"/Correlation_metadata_variables"),".pdf")
  #png(file=filename,width=600, height=550, units='mm', res=300)
  pdf(file=filename,paper="A4")
  par(mfrow=c(1,1),mar=c(1,1,1,1)) 
  nums <- unlist(lapply(metadata, is.numeric))  
  corrplot(cor(metadata[,nums],method="spearman"),method="ellipse",tl.col = "black",tl.cex=1)
  dev.off()
  #graphics.off()
  
  #check outliers:
  nums <- unlist(lapply(metadata, is.numeric))
  for (i in seq(1,length(nums)))
  {
    if (nums[[i]]==T){
      filename <- paste0(paste0(paste0(path_results,"Boxplot_outlier_check_"),names(nums)[i]),".png")
      png(file=filename,width=600, height=550, units='mm', res=300)
      par(mfrow=c(1,1),mar=c(1,1,1,1)) 
      boxplot(metadata[,names(nums)[i]], col ="orange", main = names(nums)[i])
      dev.off()
      graphics.off()
    }
  }
}

#visualize spearman correlation of logFC of data set halves for varying number of RUV factors
visualize_RUV_grid_search<- function(M_cor,path_results){
  browser()
  m<-colMeans(M_cor,na.rm=TRUE)
  e_lower <- colQuantiles(M_cor,probs=0.05,na.rm=TRUE,useNames = FALSE)
  e_upper <- colQuantiles(M_cor,probs=0.95,na.rm=TRUE,useNames = FALSE)
  n_RUVs<-colnames(M_cor)
  n_RUVs <- parse_number(n_RUVs)
  pdf(file=paste0(path_results,'/Optimal_RUV_factors.pdf'), paper="A4")
  #par(mfrow=c(1,1),mar=c(1,1,1,1)) 
  plot(n_RUVs,e_lower,type="l",ylab="Spearman correlation of logFCs",xlab="Number of RUV facors",ylim=c(-0.5,max(M_cor,na.rm=TRUE)),col='grey20',lty='dashed')
  lines(n_RUVs,e_upper,type="l",col='grey20',lty='dashed')
  #lines(n_RUVs,rep(0,length(n_RUVs)), type="l",col='black',lty='dashed')
  x = c(n_RUVs,rev(n_RUVs))
  y = as.numeric(c(e_upper,rev(e_lower)))
  polygon(x=x[!is.na(y)],y=y[!is.na(y)], col=rgb(0.05, 0.3, 0.55, 0.5), border=NA)
  lines(n_RUVs,m,type="b",pch=16,col='black',lwd=1)
  title(ylab="Spearman correlation of logFCs",xlab="Number of RUV facors")
  dev.off()
  #graphics.off()
  #ggsave(file=paste0(path_results,'/Optimal_RUV_factors.pdf'))
}

plot_PCA_of_metadata <- function(metadata_imputed_scaled,cols_numeric,path_results){
  metadata_num <- metadata_imputed_scaled[,colnames(metadata_imputed_scaled) %in% cols_numeric]
  res.pca <- princomp(metadata_num,cor = TRUE)
  library(factoextra)
  #variance
  pdf(file=paste0(path_results,'/PCA_metadata_variance_explained.pdf'), paper="A4")
  print(fviz_eig(res.pca))
  dev.off()
  #plot 1st and 2nd PC:
  rownames(res.pca$scores)<- metadata_imputed_scaled$donor_ID_python
  pdf(file=paste0(path_results,'/PCA_PC1_vs_PC2.pdf'), paper="A4")
  print(fviz_pca_ind( res.pca,
                      col.ind = "cos2", # Color by the quality of representation
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE     # Avoid text overlapping
  ))
  dev.off()
  #loadings:
  pdf(file=paste0(path_results,'/PCA_absolute contribution_PC1.pdf'))
  par(mar=c(5,10,4,2)) 
  barplot(sort(abs(res.pca$loadings[,1])),horiz=TRUE,cex.names=0.8,las=2,xlab="Absolute contribution to PC 1")
  dev.off()
  pdf(file=paste0(path_results,'/PCA_absolute contribution_PC2.pdf'))
  par(mar=c(5,10,4,2)) 
  barplot(sort(abs(res.pca$loadings[,2])),horiz=TRUE,cex.names=0.8,las=2,xlab="Absolute contribution to PC 2")
  dev.off()
}

#plot RLE = log-ratio of read count to median read count across sample
plot_RUV_RLE = function(set,DDS_RUV,path_res_RUV,n_RUV_factors){
  pdf(file=paste0(path_res_RUV,paste0(paste0('/RLE_RUV_genes_',n_RUV_factors),'.pdf')))
  sampleNames(set) <- DDS_RUV$donor_ID_python
  palette("Tableau 10")
  plotRLE(set, outline=FALSE, ylim=c(-2,2), col = factor(DDS_RUV$Group), xlab = "Samples", ylab = "Log-ratio of read count to median read count")
  legend(1, 2, legend=levels(DDS_RUV$Group), fill=palette()[1:2], cex=0.8)
  dev.off()
}

#PCA based on control genes used for RUV analysis
plot_RUV_PCA = function(set,DDS_RUV,path_res_RUV,n_RUV_factors){
  pdf(file=paste0(path_res_RUV,paste0(paste0('/PCA_RUV_genes_',n_RUV_factors),'.pdf')))
  col_vec = as.character(DDS_RUV$Group)
  col_vec[col_vec=='CTRL'] = palette()[1]
  col_vec[col_vec=='SCZ'] = palette()[2]
  plotPCA(set, col = col_vec, cex=1.2) #col=colors[DDS_RUV$Group], 
  dev.off()
}

# Bar plot of variance fractions for the first 10 genes
plot_RUV_var_fractions_per_gene = function(vp,path_res_RUV,n_RUV_factors, opt_downsampled_data, i_ds){
  if (opt_downsampled_data==TRUE){
    file_ending <- paste0(as.character(i_ds),".pdf")
  } else{
    file_ending <- ".pdf"
  }
  filename <- paste0(path_res_RUV,"/RUV_factors_percent_bars_",n_RUV_factors,file_ending)
  pdf(file=filename)
  print(plotPercentBars(vp[1:10,]))
  dev.off()
}

#Figure 1b
# violin plot of contribution of each variable to total variance
plot_RUV_var_per_factor = function(vp,path_res_RUV,n_RUV_factors,opt_downsampled_data,i_ds){
  if (opt_downsampled_data==TRUE){
    file_ending <- paste0(as.character(i_ds),".pdf")
  } else{
    file_ending <- ".pdf"
  }
  filename <- paste0(paste0(paste0(path_res_RUV,"/RUV_factors_Variance_per_factor_"),n_RUV_factors),file_ending)
  pdf(file=filename)
  print(plotVarPart(vp))
  dev.off()
}

#cumulative variance:
plot_RUV_cumulative_var = function(vp,path_res_RUV,n_RUV_factors,opt_downsampled_data,i_ds){
  if (opt_downsampled_data==TRUE){
    file_ending <- paste0(as.character(i_ds),".pdf")
  } else{
    file_ending <- ".pdf"
  }
  # mean variance per factor
  var_mean <- colSums(vp)/dim(vp)[1]
  filename <- paste0(paste0(paste0(path_res_RUV,"/RUV_factors_cumulative_variance_"),n_RUV_factors),file_ending)
  pdf(file=filename)
  plot(seq(1,length(var_mean)-1),cumsum(var_mean)[1:length(var_mean)-1]*100,type="l",ylim=c(0,100),col="black",xlab="RUV factors", ylab="Cumulative variance",xaxt = "n")
  axis(1,at=seq(1,length(var_mean)-1),labels=colnames(vp)[1:length(var_mean)-1])
  dev.off()
}


#initialize results dataframe
initialize_results_dataframe <- function(){
  R <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(R) <- c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','size_factor_estimation','shrinkage','celltype')
  return(R)  
}

#store result of all genes in dataframe:
save_results_all_genes <- function(RES,R,res,shrinkage_type,sfe_type,path_res,celltype,i){
  RES[[shrinkage_type]] <- res
  write.csv(res,paste0(path_res,paste0('df_results_',shrinkage_type,as.character(i),'.csv')),row.names = TRUE)
  if (dim(res[!is.na(res$padj) & res$padj<=alpha_val,])[1]!=0){
    R <- get_significant_results(R,res,alpha_val,sfe_type,shrinkage_type,celltype)
  }
  return(list(R=R,RES=RES))
}

#select pathways:
select_pathways <- function(c,alpha_pathways,version_nr,p_val_filter_str){
  if (nrow(c)>0){
    d <- c %>% mutate(pathID = paste0("path", row_number()))
    Pbonf_var_str <- paste0("Pbonf",str_split(as.character(alpha_pathways),"0.",2)[[1]][2])
    if (version_nr == 'v1'){
      eval(parse(text= paste0(paste0(paste0(paste0("d$",Pbonf_var_str)," <- (d$Phyper < 0."),str_split(as.character(alpha_pathways),"0.",2)[[1]][2]), "/nrow(d))")))
    } else{
      eval(parse(text= paste0(paste0(paste0(paste0("d$",Pbonf_var_str)," <- (d$",p_val_filter_str,"< 0."),str_split(as.character(alpha_pathways),"0.",2)[[1]][2]),")")))
    }
    d <- eval(parse(text=paste0(paste0("filter(d, ",Pbonf_var_str),"==TRUE)")))   # include only sig pathways
  } else{
    d <- data.frame()
  }
  return(d)
}

#=== LOAD ===
# list GMLOC: LOAD()
# get data:   a <- LOAD("fn.tsv.gz")
LOAD <- function(file,path) {
  message("Usage // LOAD() to list // a<-LOAD(fn.tsv.gz) to load")
  if ( nargs()==0 ) { system2("ls", args=path) }
  else {
    D <- fread(file = paste(path, file, sep="/"))
    return(D)
  }
}

GSA_v3 <- function(background, list_column_name, genesetsA, genesetsB, min_number_genes_in_gene_set) {
  require(tidyverse)
  require(data.table)
  require(rstatix)
  message("GSA :: hypergeometric gene set analysis by pfs 01/2023")
  message("GSA :: use clean input df (no missings, all ensgid are in geneMatrix)")
  
  #keep all the variables in background that are listed (here ensgid and what is provided as DEG_list_column_name)
  list_to_be_tested <- select(background, ensgid, {{list_column_name}}) 
  
  # sanitize input (drop missing and duplicates)
  #background <- background %>% select(ensgid, {{DEG_list_to_be_tested}})
  
  TestVar <-  list_column_name

  filter_str <- paste0(paste0(paste0(paste0("filter(is.na(ensgid)==FALSE, str_sub(ensgid, 1, 4) == \"ENSG\",is.na(",list_column_name),")==FALSE,"),list_column_name)," %in% c(TRUE, FALSE))")
  select_mutate_str <- paste0(paste0(" %>% mutate(Test = ",list_column_name),") %>% select(ensgid, Test)")
  background_clean <- eval(parse(text=paste0(paste0("background %>% ",filter_str),select_mutate_str)))
  
  background_clean <- unique(background_clean, by = "ensgid")
  background_clean_sum <- sum(background_clean$Test)
  message(paste("GSA ::", nrow(background), "rows input,", nrow(background_clean), "after cleaning"))
  message(paste("GSA :: for", TestVar, "fraction TRUE", background_clean_sum/nrow(background_clean)))
  
  # read gene sets (are unique), genesetID
  #yy <- fread(paste0(GMLOC, "/genesets.tsv"))  ### CHANGE to fit your environment
  
  # left join to preserve all genesets
  # assume :: all rows in input are valid (no missing or blank values) :: if not, nrow(background_clean) will be too large
  # note :: 1-2% of geneset ensgid will not be in input (eg input is PC, genesets have lncRNA)
  # compute worst case Phyper too
  Genesets <- merge(gene_setsB, background_clean, by.x = "ensgid", by.y = "ensgid", all.x = TRUE, all.y = FALSE)

  # counting
  Test <- Genesets %>%
    mutate(Test2 = Test) %>% 
    mutate(Test2 = replace_na(Test2,FALSE)) %>%
    #mutate(Test2 = replace_na(Test2,value=FALSE)) %>%
    #replace_na(list(Test2=FALSE)) %>% #was a necessary change since replace_na function returned an error
    group_by(group, subgroup, genesetID) %>% 
    summarise(nVar  = sum(is.na(Test)==FALSE),
              xVar  = sum(Test, na.rm = T),
              nVarX = n(),
              xVarX = sum(Test2),
              nVar_ensgid = paste(unique(ensgid[(Test==TRUE & !is.na(ensgid))]),collapse='_')) %>% #ensgids of all genes in gene set, want ensgids of the genes in test set 
    ungroup %>%
    filter(xVar >= min_number_genes_in_gene_set) %>%   # minimum size
    mutate(TestVar = TestVar,
           Ngene = nrow(background_clean),   # background is intended post-QC, not nrow(background)
           Rflag = background_clean_sum,
           Phyper = phyper(xVar-1, Rflag, Ngene-Rflag, nVar, lower.tail = FALSE),
           phCheck = phyper(xVarX-1, Rflag, nrow(background)-Rflag, nVarX, lower.tail = FALSE))
  message(paste("GSA ::", nrow(Test), "gene set tests, minimum overlap", min_number_genes_in_gene_set, "genes"))
  
  
  # add significance corrections
  if (nrow(Test)>0){
    xx <- Test %>% 
      adjust_pvalue(p.col = "Phyper", output.col = "P.fdr.all", method = "fdr") %>% 
      adjust_pvalue(p.col = "Phyper", output.col = "P.bonf.all", method = "bonferroni") %>% 
      mutate(tests.all = nrow(Test))
    ww <- xx %>% 
      group_by(group) %>% 
      summarise(tests.group = n()) %>% 
      ungroup
    vv <- xx %>% 
      group_by(group) %>% 
      adjust_pvalue(p.col = "Phyper", output.col = "P.fdr.group", method = "fdr") %>% 
      group_by(group) %>% 
      adjust_pvalue(p.col = "Phyper", output.col = "P.bonf.group", method = "bonferroni") %>% 
      ungroup %>% 
      left_join(ww) %>% 
      #rename(P.hyper = Phyper) %>% 
      mutate(P.fdr05.group = P.fdr.group < 0.05,
             P.bonf05.group = P.bonf.group < 0.05)
    Results <- genesetsA %>% 
      inner_join(vv) %>%
      select(TestVar, group, subgroup, geneset,
             genes.in.backround=Ngene, genes.TestVar.true=Rflag, genes.in.geneset=nVar, overlap.TestVar.geneset=xVar, genes.in.geneset.tested=nVar_ensgid,
             Phyper, 
             ends_with(".group"), ends_with(".all"),
             phCheck, xVarX, nVarX) %>% 
      arrange(Phyper, by.group=TRUE) 
    rm(background_clean,background_clean_sum,genesetsA,genesetsB,xx,ww,TestVar,Genesets,Test)
  } else{
    Results <- data.frame()
    rm(background_clean,background_clean_sum,genesetsA,genesetsB,TestVar,Genesets,Test)
  }
  
  return(Results)
  
}

save_pathways_as_csv_file<- function(result_GSA_pathways,result_GSA_sign_pathways, path_GSA_results,CT_str,info_str){
  #if there are pathways found, save them in excel sheets:
  if (dim(result_GSA_pathways)[1]>0){
    #make sure results path exists:
    if (dir.exists(path_GSA_results)==FALSE){
      dir.create(path_GSA_results,recursive=TRUE,showWarnings = FALSE)
    }
    setwd(path_GSA_results)
    
    file_ending_all <- paste0(info_str,"_all.csv")
    if (dim(result_GSA_sign_pathways)[1]>0){
      file_ending_sign <- paste0(info_str,"_sign.csv")
    }
    
    if (!is.null(CT_str)){
      write.csv2(result_GSA_pathways,paste0(CT_str,"_",file_ending_all),row.names=FALSE)
      if (dim(result_GSA_sign_pathways)[1]>0){
        write.csv2(result_GSA_sign_pathways,paste0(CT_str,"_",file_ending_sign),row.names=FALSE)
      }
    }
  }
}

run_pathway_clustering_with_rrvgo <- function(opt_input_data, settings, main_path,p_val_cutoff){
  if (str_detect(opt_input_data,"GRN")){
    settings$n_clusters_range <- get_cell_type_clusters_GRN(settings$opt_metacell_settings)
  }
  for (n_clusters in settings$n_clusters_range){
    if (opt_input_data == "DEG_GSA_result" || opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result"){  
      if (n_clusters==3){
        cluster_folder <- paste0(n_clusters,'_classes/')
      } else{
        cluster_folder <- paste0(n_clusters,'_CTs/')
      }
    } else if (opt_input_data=="GRN_modules_gprofiler_result"){
      cluster_folder <- paste0(n_clusters,'/')
    } else{
      cluster_folder <- ""
    }
    for (sub_folder in settings$sub_folders){
      for (min_overlap in settings$min_overlap_range){
        for (background_folder in settings$background_folders){
          
          #load significant pathways:
          data_and_results_path <- get_results_path(main_path,opt_input_data,cluster_folder,min_overlap,background_folder,sub_folder,settings$opt_metacell_settings)
          
          D <- get_data(opt_input_data,main_path, data_and_results_path,p_val_cutoff, cf, sf, settings$opt_up_and_down_sep,settings$opt_metacell_settings)
          
          if (settings$opt_up_and_down_sep==FALSE || opt_input_data=="DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result"){
            crit <- is.null(D)==FALSE && (dim(D)[1]>0)
            reg_strings <- c('all')
          } else{
            crit <- is.null(D$up)==FALSE && (dim(D$up)[1]>0) && is.null(D$down)==FALSE && (dim(D$down)[1]>0)
            reg_strings <- c('up','down')
          }
          if (crit){
            #extract group, subgroup, geneset, 1-P.bonf.group as score
            #put together all up-reg pathways across CTs and all down-reg pathways across CTs 
            for (reg_str in reg_strings){
              for (GO_category in c('BP','MF','CC')){
                #filter data set for relevant rows (pathways)
                if (settings$opt_up_and_down_sep==FALSE || opt_input_data=="DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result"){
                  D_sel <- D %>% filter(subgroup==paste0("GO:",GO_category))
                } else {
                  if (reg_str == 'up'){
                    D_sel <- D$up %>% filter(subgroup==paste0("GO:",GO_category))
                  } else {
                    D_sel <- D$down %>% filter(subgroup==paste0("GO:",GO_category))
                  }
                } 
                if (all(colnames(D_sel)!="ID")){
                  D_sel$ID <- D_sel %>% {str_extract(.$geneset, "^.{10}")}
                }
                
                #update settings$CT_groups to get all modules in case of separate module analysis of GRN result:
                if (opt_input_data == "GRN_modules_gprofiler_result"){
                  if (opt_modules_separately ==TRUE){
                    settings$CT_groups <- unique(D_sel$Module)
                  }
                }
                for (ct_group in settings$CT_groups){
                  if (opt_input_data == "GRN_modules_gprofiler_result"){
                    if (opt_modules_separately ==TRUE){
                      D_sel_ct <- D_sel %>% filter(grepl(ct_group, Module))
                    } else{
                      D_sel_ct <- D_sel # nothing filtered, put all modules together
                    }
                  } else{
                    #select current cell type (group)
                    if (ct_group=='all'){
                      D_sel_ct <- D_sel # nothing filtered
                    } else if (ct_group=='Neurons'){
                      D_sel_ct <- D_sel %>% filter(grepl('*neurons*',TestVar))
                    } else{
                      D_sel_ct <- D_sel %>% filter(str_detect(TestVar,paste0('^',ct_group,'_')))
                    }
                  }
                  
                  
                  if (dim(D_sel_ct)[1]>2){
                    #load_OrgDb("org.Hs.eg.db")
                    
                    # 1) calculate score similarity matrix
                    simMatrix <- calculateSimMatrix(D_sel_ct$ID,orgdb="org.Hs.eg.db",ont=GO_category,method="Rel")
                    
                    if (!is.null(dim(simMatrix))){
                      if (settings$opt_sim_based_on == 'size'){
                        reducedTerms <- reduceSimMatrix(simMatrix,
                                                        threshold=0.7,
                                                        orgdb="org.Hs.eg.db")
                      } else if (settings$opt_sim_based_on == 'scores'){
                        scores <- setNames(-log10(D_sel_ct$P.bonf.group), D_sel_ct$ID)
                        reducedTerms <- reduceSimMatrix(simMatrix,
                                                        scores,
                                                        threshold=0.7,
                                                        orgdb="org.Hs.eg.db")
                      }
                      #save table of reduced terms:
                      if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
                        filename_RT <- paste0('reduced_terms_rrvgo_',ct_group,'_',GO_category,'_',settings$opt_sim_based_on,'_',reg_str,'DEG_GSA_and_GRN_modules_gprofiler.csv')
                        write.csv(reducedTerms,paste0(data_and_results_path,filename_RT))
                      } else if (opt_input_data == "all_GRN_modules_gprofiler_result"){
                        filename_RT <- paste0('reduced_terms_rrvgo_',ct_group,'_',GO_category,'_',settings$opt_sim_based_on,'_',reg_str,'_GRN_modules_gprofiler.csv')
                        write.csv(reducedTerms,paste0(data_and_results_path,filename_RT))
                      } else{
                        filename_RT <- paste0('reduced_terms_rrvgo_',ct_group,'_',GO_category,'_',settings$opt_sim_based_on,'_',reg_str,'.csv')
                        write.csv(reducedTerms,paste0(data_and_results_path,filename_RT))
                        
                        module_list <- plot_pathway_clustering(simMatrix, reducedTerms, data_and_results_path, settings$opt_sim_based_on, GO_category, reg_str, ct_group, settings$opt_up_and_down_sep, opt_input_data)
                        filename_ML <- paste0('module_list_rrvgo_',ct_group,'_',GO_category,'_',settings$opt_sim_based_on,'_',reg_str,'.csv')
                        write.csv(reducedTerms,paste0(data_and_results_path,filename_ML))
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

save_GO_ensgids <- function(df1,search_terms,search_term_IDs,GO_str){
  for (term_i in seq(1,length(search_term_IDs))){
    print(search_terms[term_i])
    i <- which(df1$term==search_term_IDs[term_i])
    tmp.genelist <- df1[i,"gene"] %>% 
      as.data.frame() 
    colnames(tmp.genelist)[1]="ENSGID"
    #export tmp.genelist:
    fwrite(tmp.genelist,file=paste0(path.out,"ensgid_GO_",GO_str,"_",search_terms[term_i],".tsv"), sep = "\t", col.names = T)
  }
}