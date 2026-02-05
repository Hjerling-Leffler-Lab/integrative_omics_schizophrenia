# Differential gene expression analysis with DESeq2
# performs likelihood ratio test for models with RUV variables
# full model: 6 RUV factors and group status
# reduced model: 6 RUV factors 
#
# date: 10.02.2021; last update 12.05.2025
# author: Lisa Bast
# version: 0.1.6

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

opt_server <- FALSE
if (opt_server ==TRUE){
  n_cores <- 12 #24
} else{
  n_cores <- 4 #12 #24
}

opt_pagoda <- TRUE#FALSE#
if (opt_pagoda==TRUE){
  add_str_pagoda <- '_pagoda'
} else{
  add_str_pagoda <- ''
}

n_RUV_vars_to_inlcude <- 6

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
path_code = getwd()
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

opt_celltype_groups = 'conos_cluster_based'

#should shrinkage be performed?:
opt_perform_shrinkage = FALSE #TRUE

#paths for data and results:
data_folder <- 'CT_clustered_aggregated_data/'

opt_aggregation = 'sum'
opt_export_sfe_normalized_data=FALSE#TRUE

design_folder_name <- paste0('design_with_',n_RUV_vars_to_inlcude,'RUV_covariates_',opt_aggregation)

filename_data_sample_info <- "T1_basic_donor_information.xlsx"
filename_data_ancestry <- "sample_brn_ancestry.tsv"

paths <- c()
paths$data <- paste0(main_path,main_subfolder,"output/",data_folder)
paths$result <- paste0(main_path,main_subfolder,"/output/DEGs/")
paths$res_RUV <- paste0(main_path,main_subfolder,"/output/RUV/",design_folder_name,"/")

#create folder if not existent and change file permissions
create_folder_and_set_permissions(paths$data)
create_folder_and_set_permissions(paths$result)
create_folder_and_set_permissions(paths$res_RUV)

#read data
setwd(paths$data)

countData_filenames <- Sys.glob(paste0("Agg_counts_",add_str_pagoda,"_TH_and_D_adj_filtered",'_*_',opt_aggregation,'.csv'))
colData_filenames <- Sys.glob(paste0("G_info_",add_str_pagoda,"_TH_and_D_adj_filtered",'_*_',opt_aggregation,'.csv'))

paths$xy_chr_list <- paste0(main_path,'/4_data_integration_and_cell_type_annotation/output/')
paths$sample_info_data <- paste0(main_path,"/2_alignment/output/")

#initialize Results data frames:
R <- initialize_results_dataframe()

for (file_id in seq(1,length(countData_filenames))){
  #load data:
  setwd(paths$data)
  data <- load_count_data(countData_filenames,colData_filenames,file_id,opt_aggregation,paths)
  data$metadata <- get_metadata(data, paths,filename_data_sample_info,filename_data_ancestry, filename_QC_data, TRUE, TRUE, TRUE)
  
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
  
  #specify results_folder
  path_res <- get_results_folder(paths$result,design_folder_name,data$results_subfolder,"",FALSE)
  create_folder_and_set_permissions(path_res)
  
  #remove genes with low counts and located on x and y chromosomes and store data in DDS object:
  DDS <- get_filtered_data(data$CC,data$metadata,design_formula,paths$xy_chr_list,p_samples_larger_10_counts)
  
  
  #determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
  sfe_type <- get_sfe_type(DDS)
  
  #run DSEq2 for very basic model withour any RUV covariates:
  DDS <- DESeq(object=DDS,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = ~ 1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
  

  # plot size factors (distribution for each group)--> check if roughly the same across groups
  plot_size_factors(DDS,data$celltype,path_res,sfe_type)

  
  #obtain results:
  #Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning data-driven weights to each hypothesis.
  #res <- results(DDS,alpha=alpha_val,filterFun=ihw,cooksCutoff=FALSE)
  res <- results(DDS,alpha=alpha_val,contrast=c("Group","SCZ", "CTRL"), pAdjustMethod = "BH", independentFiltering = FALSE)
  #summary(res)
  #how many adjusted p-values were less than alpha_val?
  #sum(res$padj < alpha_val, na.rm=TRUE)

  #store result of all genes in dataframe results_obj$RES and result of significant genes in seperate dataframe (results_obj$R)
  RES = list()
  
  results_obj <- get_results_object_all_genes(RES,R,res,shrinkage_types[1],sfe_type,path_res,data$celltype,c(),alpha_val)
  
  plot_significant_genes(DDS,res,alpha=alpha_val,mode=shrinkage_types[1],path_res)
  
  #save dataframe with all significant genes
  if (dim(results_obj$R)[1]>0){
	  write.csv(results_obj$R,paste0(path_res,'df_results_',opt_aggregation,'.csv'),row.names = TRUE)
  }
  
  #shrinkage: yes/ no?
  if (opt_perform_shrinkage){
		for (st in shrinkage_types){
		  if (st != "no_shrinkage"){
  			res <- lfcShrink(DDS, coef="Group_SCZ_vs_CTRL", type=st)
  			#summary(res)
  			results_obj <- get_results_object_all_genes(results_obj$RES,results_obj$R,res,st,sfe_type,path_res,data$celltype,c(),alpha_val)
  			plot_significant_genes(DDS,res,alpha=alpha_val,mode=st,path_res)
		  }
		}
		#plot MA plots (all shrinkage types)
		plot_MA_for_various_shrinkage_types(results_obj$RES,shrinkage_types,path_res)
  }
  

  if ((n_RUV_vars_to_inlcude>0)&(n_RUV_vars_to_inlcude<=8)){
	  print(paste0("n_RUV_vars_to_inlcude: ",n_RUV_vars_to_inlcude))
    design_folder_name <- paste0('design_with_',n_RUV_vars_to_inlcude,'_RUVs')
    #update results_folder
    path_res <- get_results_folder(paths$result,design_folder_name,data$results_subfolder,"",FALSE)
    create_folder_and_set_permissions(path_res)
		#perform RUV analysis to test of design should be more complicated or if simple design is sufficient
		if (n_RUV_vars_to_inlcude>0){
			#load common factors:
			setwd(paths$res_RUV)
      filename_RUV_factors <- paste0("RUV_factors_",n_RUV_vars_to_inlcude,"_factors.Rds")
      if (!(file.exists(filename_RUV_factors))){
			  print(paste0("file ",filename_RUV_factors," does not exist in directory ",paths$res_RUV))
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
			  design_str <- get_design_str(factors_to_include,TRUE,FALSE)
			  design_str_red <- get_design_str(factors_to_include,TRUE,TRUE)
			  
			  
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
			  #update design in DDS object:
			  DDS_RUV@design <- eval(parse(text=design_str))
			  DDS_RUV <- DESeq(object=DDS_RUV,parallel=opt_parallel,sfType = sfe_type, test="LRT", reduced = eval(parse(text=design_str_red)), useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
			  #save results:
			  res_RUV <- results(DDS_RUV,alpha=alpha_val,contrast=c("Group", "SCZ", "CTRL"), pAdjustMethod = "BH", independentFiltering = FALSE)
			  RES_RUV = list()
			  
			  results_obj_RUV <- get_results_object_all_genes(RES_RUV,R,res_RUV,shrinkage_types[1],sfe_type,path_res,data$celltype,c(),alpha_val)
			  
			  if (n_RUV_vars_to_inlcude>0 & dim(results_obj_RUV$R)[1]>0){
				  write.csv(results_obj_RUV$R,paste0(path_res,'df_results_',opt_aggregation,'_RUV','.csv'),row.names = TRUE)
			  }

			  #shrinkage: yes/ no?
			  if (opt_perform_shrinkage){
					for (st in shrinkage_types){
					  if (st != "no_shrinkage"){
							res_RUV <- lfcShrink(DDS_RUV, coef="Group_SCZ_vs_CTRL", type=st)
							#summary(res)
							results_obj_RUV <- get_results_object_all_genes(results_obj_RUV$RES,results_obj_RUV$R,res_RUV,st,sfe_type,path_res,data$celltype,c(),alpha_val)
							plot_significant_genes(DDS_RUV,res_RUV,alpha=alpha_val,mode=st,path_res)
					  }
					}
					#plot MA plots (all shrinkage types)
					plot_MA_for_various_shrinkage_types(results_obj_RUV$RES,shrinkage_types,path_res)
			  } 
			}
			
		}
  }

}

