# about: Differential gene expression analysis with DESeq2; option to run on randomized group labels
#        LR-Test, correcting for Sex, PMI, Age, Library, source (of tissue, 3 sources in total), ancestry_european 
#        size factor estimation with scran:computeSumFactors
# date: 10.02.2021; last update 30.01.2025
# author: Lisa Bast
# version: 0.2.0

library("DESeq2")
library("glmGamPoi")
library("iotools") # requires for read.csv.raw efficiently read csv files
library("BiocParallel")
library("IHW")
library("stringr")
library("ashr")
library("scuttle")
library("dplyr")
library("loomR")
library("scran")
library("readxl")
library("mice")

set.seed(17274)
graphics.off()

opt_server<-FALSE

n_cores <- 4 #12 #24

opt_pagoda <- TRUE#FALSE#
if (opt_pagoda==TRUE){
  add_str_pagoda <- '_pagoda'
} else{
  add_str_pagoda <- ''
}

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
main_subfolder = "/5b_differential_gene_expression_analysis/"


alpha_val <- 0.05
opt_filter_genes_with_low_counts <- TRUE
if (opt_filter_genes_with_low_counts){
  p_samples_larger_10_counts<-70
}

opt_ancestry_binary <- TRUE
opt_scale_and_center_metadata <- TRUE
opt_load_metadata <- TRUE # set to true once they are put together

#opt_celltype_groups = 'scmap'
opt_celltype_groups = 'conos_cluster_based'

#should shrinkage be performed?:
opt_perform_shrinkage = FALSE #TRUE

#paths for data and results:
data_folder <- 'CT_clustered_aggregated_data/'

opt_aggregation = 'mean'

filename_data_sample_info <- "T1_basic_donor_information.xlsx"
filename_data_ancestry <- "sample_brn_ancestry.tsv"
filename_QC_data <- "T2_quality_metrics_per_sample.xlsx"
filename_scdata <- "Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered_subsampled_to_10_percent_cells.loom"

paths <- get_paths(main_path,main_subfolder,data_folder)

#read data

setwd(paths$data)
countData_filenames <- Sys.glob(paste0("Agg_counts_",add_str_pagoda,"_TH_and_D_adj_filtered_",'*_',opt_aggregation,'.csv'))
colData_filenames <- Sys.glob(paste0("G_info_",add_str_pagoda,"_TH_and_D_adj_filtered_",'*_',opt_aggregation,'.csv'))

#initialize Results data frames:

R <- initialize_results_dataframe()

for (file_id in seq(1,length(countData_filenames))){
  #load data:
  setwd(paths$data)
  data <- load_count_data(countData_filenames,colData_filenames,file_id,opt_celltype_groups,opt_aggregation,paths)
  data$metadata <- get_metadata(data, paths,filename_data_sample_info,filename_data_ancestry, filename_QC_data, opt_load_metadata, opt_scale_and_center_metadata, opt_ancestry_binary)
  
  #specify design:
  design_formula <- formula(~ Sex + PMI_h + Age + Library+ source + ancestry_european + Group)
  design_folder_name <- c('design_with_6_covariates')
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
  path_res <- get_results_folder(paths$results,design_folder_name,data$results_subfolder,"",FALSE)
  
  #remove genes with low counts and located on x and y chromosomes and store data in DDS object:
  DDS <- get_filtered_data(data$CC,data$metadata,design_formula,paths$xy_chr_list)
  
  #scran::computeSumFactors
  #load single cell data to calculate sizeFactors with scran
  counts <- get_non_aggregated_counts(paths$scdata,filename_scdata)
  browser()
  sizeFactors(DDS) <- scran::computeSumFactors(x=counts)
  
  #determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
  sfe_type <- get_sfe_type(DDS)
  
  #run DSEq2:
  DDS <- DESeq(object=DDS,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = ~ 1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
  
  
  
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
  



