# Differential gene expression analysis; including sophisticated models with RUV and MARS variables; option to run for randomized group labels
# with DESeq2
# date: 10.02.2021; last update 21.03.2022
# author: Lisa Bast
# version: 0.1.5

library("DESeq2")
library("glmGamPoi")
library("iotools") # requires for read.csv.raw efficiently read csv files
library("BiocParallel")
library("IHW")
library("stringr")
library("ashr")
library("scuttle")
library("dplyr")
library("RUVSeq")#RUV analysis
library("loomR")

set.seed(17274)
graphics.off()

n_cores <- 4 #12 #24

n_cluster_range <- c(16) #c(16,37)

opt_pagoda <- TRUE#FALSE#
if (opt_pagoda==TRUE){
  add_str_pagoda <- '_pagoda'
} else{
  add_str_pagoda <- ''
}

n_RUV_vars_to_inlcude_range <- c(0,4,6,8) 

#settings for parallelization:
if (n_cores>1){
  opt_parallel = TRUE
  if (opt_server==TRUE){
    register(MulticoreParam(n_cores))
  } else{
    register(SnowParam(n_cores))
  }
} else{
  opt_parallel=FALSE
}

#define pathes and load functions
code_path = getwd()
source("utils.R")

setwd("../")
setwd("../")
main_path = getwd()
main_subfolder = "5b_differential_gene_expression_analysis/"


alpha_val <- 0.05
opt_filter_genes_with_low_counts <- TRUE
if (opt_filter_genes_with_low_counts){
  p_samples_larger_10_counts<-70
}

#opt_celltype_groups = 'scmap'
opt_celltype_groups = 'conos_cluster_based'

#should shrinkage be performed?:
opt_perform_shrinkage = FALSE #TRUE

#paths for data and results:
data_folder <- 'CT_clustered_aggregated_data/'

opt_aggregation = 'sum'
opt_export_sfe_normalized_data=FALSE#TRUE

filename_data_sample_info <- "T1_basic_donor_information.xlsx"
filename_data_ancestry <- "sample_brn_ancestry.tsv"

for (n_cluster in n_cluster_range){
	if (n_cluster==3){
		opt_calculate_DEGs_for = 'classes'
	} else{
		opt_calculate_DEGs_for = 'cell types'
	}
  
	if (opt_calculate_DEGs_for == 'cell types'){
	  path_data <- paste0(main_path,main_subfolder,"/data/",data_folder,"/",n_cluster,"_CTs/")
	  path_results <- paste0(main_path,main_subfolder,"/output/DEGs/",n_cluster,'_CTs/')
	} else{
	  path_data <- paste0(main_path,main_subfolder,"/data/",data_folder,"/",'3_classes/')
	  path_results <- paste0(main_path,main_subfolder,"/output/DEGs/3_classes/")
	}
	path_res_RUV <- paste0(main_path,main_subfolder,"/output/RUV/")
	
	#read data
	setwd(path_data)

  countData_filenames <- Sys.glob(paste0("Aggregated_counts",add_str_align,add_str_pagoda,"_TH_and_D_adj_filtered_",opt_celltype_groups,'_*_',opt_aggregation,'.csv'))
  colData_filenames <- Sys.glob(paste0("G_info",add_str_align,add_str_pagoda,"_TH_and_D_adj_filtered_",opt_celltype_groups,'_*_',opt_aggregation,'.csv'))

	path_xy_chr_list <- paste0(main_path,'4_data_integration_and_cell_type_annotation/output/')
	path_sample_info_data <- paste0(main_path,"2_alignment/output/")
	
	#initialize Results data frames:
	R <- initialize_results_dataframe()

	ds_file_ids <- sapply(strsplit(sapply(strsplit(countData_filenames,'_ds_'),"[",2),".csv"),"[",1) 
	for (file_id in seq(1,length(countData_filenames))){
	  #load data:
	  setwd(path_data)
	  data <- load_count_data(countData_filenames,colData_filenames,file_id,opt_celltype_groups,opt_aggregation,path_data, add_str_align)
	  data$metadata <- get_metadata(data, path_sample_info_data,filename_data_sample_info,filename_data_ancestry, path_QC_data, filename_QC_data, c(''), TRUE, TRUE, add_str_align)

	  #specify design:
	  design_formula <- formula(~ Group)
	  design_folder_name <- c('design_without_covariates')
	  opt_design_with_interactions <- FALSE
	  
	  #specify shrinkage of effect size (log fold change estimates):
	  if (opt_perform_shrinkage){
  		if (opt_design_with_interactions==TRUE){
  		  shrinkage_types=c("no_shrinkage", "apeglm","ashr")
  		} else{
  		  shrinkage_types=c("no_shrinkage", "apeglm","ashr","normal")
  		}
	  }else{
		  shrinkage_types=c("no_shrinkage")
	  }
	  
	  #secify results_folder
	  path_res <- get_results_folder(path_results,design_folder_name,data$results_subfolder,opt_use_downsampled_data,"",FALSE)
	  
	  #remove genes with low counts and located on x and y chromosomes and store data in DDS object:
	  DDS <- get_filtered_data(data$CC,data$metadata,design_formula,path_xy_chr_list)
	  
	  #determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
	  sfe_type <- get_sfe_type(DDS)
	  
	  #run DSEq2:
	  DDS <- DESeq(object=DDS,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = ~ 1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
	  
	  #export filtered and sfe normalized data for wilcoxon to .loom file
	  export_sfe_normalized_data(opt_export_sfe_normalized_data,DDS,path_data,countData_filenames,file_id)
	  
	  #write.csv(data_sfe_normalized,new_filename, row.names = TRUE, header = paste(colData(DDS)$donor_ID_python,colData(DDS)$Group,sep='_'))
	  
	  # plot size factors (distribution for each group)--> check if roughly the same across groups
	  plot_size_factors(DDS,data$celltype,path_res,sfe_type,i)

	  
	  #obtain results:
	  #Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning data-driven weights to each hypothesis.
	  #res <- results(DDS,alpha=alpha_val,filterFun=ihw,cooksCutoff=FALSE)
	  res <- results(DDS,alpha=alpha_val,contrast=c("Group","SCZ", "CTRL"), pAdjustMethod = "BH", independentFiltering = FALSE)
	  #summary(res)
	  #how many adjusted p-values were less than alpha_val?
	  #sum(res$padj < alpha_val, na.rm=TRUE)
	  
	  #store result of all genes in dataframe results_obj$RES and result of significant genes in seperate dataframe (results_obj$R)
	  RES = list()
	  results_obj <- save_results_all_genes(RES,R,res,shrinkage_types[1],sfe_type,path_res,data$celltype,i)
	  
	  plot_significant_genes(DDS,res,alpha=alpha_val,mode=shrinkage_types[1],path_res,i)
	  
	  #save dataframe with all significant genes
	  if (dim(results_obj$R)[1]>0){
		  write.csv(results_obj$R,paste0(path_res,'df_results_',opt_aggregation,as.character(i),'.csv'),row.names = TRUE)
	  }
	  
	  #shrinkage: yes/ no?
	  if (opt_perform_shrinkage){
  		for (st in shrinkage_types){
  		  if (st != "no_shrinkage"){
    			res <- lfcShrink(DDS, coef="Group_SCZ_vs_CTRL", type=st)
    			#summary(res)
    			results_obj <- save_results_all_genes(results_obj$RES,results_obj$R,res,st,sfe_type,path_res,data$celltype,i)
    			plot_significant_genes(DDS,res,alpha=alpha_val,mode=st,path_res)
  		  }
  		}
  		#plot MA plots (all shrinkage types)
  		plot_MA_for_various_shrinkage_types(results_obj$RES,shrinkage_types,path_res)
	  }
	  
	  for (n_RUV_vars_to_inlcude in n_RUV_vars_to_inlcude_range){
		    if ((n_RUV_vars_to_inlcude>0)&(n_RUV_vars_to_inlcude<=8)){
		      if (n_RUV_vars_to_inlcude>0){
		        design_folder_name <- paste0('design_with_',n_RUV_vars_to_inlcude,'_RUVs')
		      }
  		    #create and get results folder:
  		    path_res_0RUV <-get_results_folder(path_results,design_folder_name,data$results_subfolder,opt_use_downsampled_data,"",FALSE)
		    }
		    if ((n_RUV_vars_to_inlcude>0)&(n_RUV_vars_to_inlcude<=8)){
  			  print(paste0("n_RUV_vars_to_inlcude: ",n_RUV_vars_to_inlcude))
  				#perform RUV analysis to test of design should be more complicated or if simple design is sufficient
  				if (n_RUV_vars_to_inlcude>0){
  					#load common factors:
  					setwd(path_res_RUV)
            filename_RUV_factors <- paste0("RUV_factors_",n_RUV_vars_to_inlcude,"_factors.Rds")
            if (!(file.exists(filename_RUV_factors))){
  					  print(paste0("file ",filename_RUV_factors," does not exist in directory ",path_res_RUV))
  					} else{
  					  RUV_factor_info <- readRDS(file=filename_RUV_factors)
  					  if (n_RUV_vars_to_inlcude == 4){
  						  factors_to_include <- c("W1","W2","W3","W4")
  					  } else if (n_RUV_vars_to_inlcude == 6){
  						  factors_to_include <- c("W1","W2","W3","W4","W5","W6")
  					  } else if (n_RUV_vars_to_inlcude == 8){
  						  factors_to_include <- c("W1","W2","W3","W4","W5","W6","W7","W8")
  					  } else{
  						  print("This option is not implemented (yet)!")
  					  }
  					  
  					  #change DDS object for current cell type
  					  DDS_RUV <- DDS
  					  s_id_Wi<-1
  					  s_id_DDS<-1
  					  sample_sel<-c()
  					  #browser()
  					  if (dim(RUV_factor_info)[1]>dim(DDS)[2]){
    						#not all samples contribute to current cell type (cell type is not present in one or some samples)
    						for (sample_Wi in rownames(RUV_factor_info)){
    						  if (sample_Wi == DDS_RUV$donor_ID_python[s_id_DDS]){
      							sample_sel[s_id_DDS] <- s_id_Wi
      							s_id_DDS=s_id_DDS+1
    						  }
    						  s_id_Wi=s_id_Wi+1
    						}
  					  } else{
  						  sample_sel <- seq(1,dim(DDS)[2])
  					  }
  					  
  					  
  					  #specify design:
  					  design_str <- get_design_str(factors_to_include,sample_sel,TRUE,FALSE)
  					  design_str_red <- get_design_str(factors_to_include,sample_sel,TRUE,TRUE)
  					  
  					  # add Wi to colData in DDS_RUV object:
  					  #browser()
  					  RUV_factor_info_selection <- RUV_factor_info[rownames(RUV_factor_info) %in% DDS_RUV$donor_ID_python,]
  					  DDS_RUV <- integrate_RUV_factors_in_DDS_obj(DDS_RUV,RUV_factor_info_selection,'','none') 
  					  
  					  eval(parse(text=design_str))
  					  
  					  print("celltype: ")
  					  cat(data$celltype)
  					  print(paste("file_id: ",file_id))
  					  
  					  opt_design_with_interactions <- FALSE
  					  #run DSEq2 with adapted design formula:
  					  #browser()
  					  if (length(DDS_RUV@colData@listData$donor_ID_python)!=length(DDS_RUV@colData@listData$W3)){
  					    next
  					  }
  					  DDS_RUV <- DESeq(object=DDS_RUV,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = eval(parse(text=design_str_red)), useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
  					  #save results:
  					  res_RUV <- results(DDS_RUV,alpha=alpha_val,contrast=c("Group", "SCZ", "CTRL"), pAdjustMethod = "BH", independentFiltering = FALSE)
  					  RES_RUV = list()
  					  results_obj_RUV <- save_results_all_genes(RES_RUV,R,res_RUV,shrinkage_types[1],sfe_type,path_res_RUV,data$celltype,i)
  					  
  					  if (n_RUV_vars_to_inlcude>0 & dim(results_obj_RUV$R)[1]>0){
  						  write.csv(results_obj_RUV$R,paste0(path_res_RUV,'df_results_',opt_aggregation,'_RUV',as.character(i),'.csv'),row.names = TRUE)
  					  }
  					  
  					  #shrinkage: yes/ no?
  					  if (opt_perform_shrinkage){
    						for (st in shrinkage_types){
    						  if (st != "no_shrinkage"){
      							res_RUV <- lfcShrink(DDS_RUV, coef="Group_SCZ_vs_CTRL", type=st)
      							#summary(res)
      							results_obj_RUV <- save_results_all_genes(results_obj_RUV$RES,results_obj_RUV$R,res_RUV,st,sfe_type,path_res_RUV,data$celltype,i)
      							plot_significant_genes(DDS_RUV,res_RUV,alpha=alpha_val,mode=st,path_res_RUV)
    						  }
    						}
    						#plot MA plots (all shrinkage types)
    						plot_MA_for_various_shrinkage_types(results_obj_RUV$RES,shrinkage_types,path_res_RUV)
  					  }
  					}
					
  				#no RUV analysis:
  				} else { 
  				  sample_sel <- seq(1,dim(DDS)[2])
  					#specify design:
  					design_str <- get_design_str(factors_to_include,sample_sel,FALSE,FALSE)
  					design_str_red <- get_design_str(factors_to_include,sample_sel,FALSE,TRUE)
  					
  					eval(parse(text=design_str))
  					
  					opt_design_with_interactions <- FALSE
  					
  					#run DSEq2 with adapted design formula:
  					DDS <- DESeq(object=DDS,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = eval(parse(text=design_str_red)), useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
  					#save results:
  					res_0RUV <- results(DDS,alpha=alpha_val,contrast=c("Group", "SCZ", "CTRL"), pAdjustMethod = "BH", independentFiltering = FALSE)
  					RES_0RUV = list()
  					results_obj_0RUV <- save_results_all_genes(RES_0RUV,R,res_0RUV,shrinkage_types[1],sfe_type,path_res_MARS,data$celltype,"")
  					
  					if (dim(results_obj_0RUV$R)[1]>0){
  					  write.csv(results_obj_0RUV$R,paste0(path_res_0RUV,'df_results_',opt_aggregation,'_0RUV','.csv'),row.names = TRUE)
  					}
  					
  					#shrinkage: yes/ no?
  					if (opt_perform_shrinkage){
  					  for (st in shrinkage_types){
    						if (st != "no_shrinkage"){
    						  res_0RUV <- lfcShrink(DDS, coef="Group_SCZ_vs_CTRL", type=st)
    						  #summary(res)
    						  results_obj_0RUV <- save_results_all_genes(results_obj_0RUV$RES,results_obj_0RUV$R,res_MARS,st,sfe_type,path_res_0RUV,data$celltype,i)
    						  plot_significant_genes(DDS,res_0RUV,alpha=alpha_val,mode=st,path_res_0RUV)
    						}
  					  }
  					  #plot MA plots (all shrinkage types)
  					  plot_MA_for_various_shrinkage_types(results_obj_0RUV$RES,shrinkage_types,path_res_0RUV)
  					}
  					
  				}
  		  }
		  }
	  }
  }



