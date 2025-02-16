#about: Build DESeq2 model by using RUV variables
#       prior to running DEG_analysis_DESeq2.R
#       Fig. S5 A-F
#author: Lisa Bast
#date: 30.11.2021
#version: 0.0.4

library('DESeq2')
library("RUVSeq")#RUV analysis
library('earth')
library('readxl')
library('VIM')
library('mice')
library("corrplot")
library("caret")
library('sjmisc')
library("tidyverse")#for transforming matrix into tibble
library("foreach")#parallel for loops
library("doParallel")
library("iotools") # requires for read.csv.raw efficiently read csv files

graphics.off()

opt_server <- TRUE

n_workers <- 20 #35
n_cluster <- 15 # cell type resolution

opt_pagoda <- TRUE#FALSE#
if (opt_pagoda==TRUE){
  add_str_pagoda <- '_pagoda'
} else{
  add_str_pagoda <- ''
}

opt_load_metadata <- TRUE #FALSE#  set to false in first run

if (opt_server==TRUE){
  args <- commandArgs(trailingOnly=TRUE)
}else{
  args <- "1"
}
ct_id <- as.integer(args)

#RUV settings:
opt_filter_high_CT_specificity_genes <- TRUE
if (ct_id==1){
  opt_run_RUV <- TRUE
} else{
  opt_run_RUV <- FALSE#TRUE#
}

opt_run_split_halves <- TRUE
alpha_val <- 0.5 # used to determine which genes are kept for RUV analysis; initial DESeq2 run needs to result in p-val larger than alpha value for a gene to be kept in RUV factor calculation
specificity_TH <- 0.5# 0.8
min_p_donors_groups <- 0.7 # for 40 donors in group SCZ at least 14 have to be in each dataset half; the higher the more balanced groups have to be
n_data_splits <- 1000
n_RUV_factors <- 20
#DESEQ2 settings:
alpha_val <- 0.05
opt_filter_genes_with_low_counts <- TRUE
if (opt_filter_genes_with_low_counts){
  p_samples_larger_10_counts<-70
}

#which cell type annotation files should be used?:
opt_aggregation = 'sum'
opt_ancestry_binary = TRUE
opt_parallel_for_loop = TRUE # for RUV grid search 
opt_DESeq2_parallel = FALSE # causes some errors I don't understand. better set to FALSE
opt_scale_and_center_metadata = TRUE#FALSE#
#which metadata files should be used?:
filename_data_sample_info <- "sample_info_SUN_SCZ.xlsx"
filename_QC_data <- "metrics_summary_tidy.csv" 
filename_data_ancestry <- 'sample_brn_ancestry.tsv'
filename_CT_spec <- 'speMx_Expr_15CTs.Rdata'

#define paths and load functions:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()

data_folder <- "aggregated_data/"
paths <- get_paths(main_path,main_project_path,data_folder)

if (opt_run_RUV == TRUE){
  ## 0) load data aggregated across all cell types:
  # expression data (aggregated over all cells from a specific sample for each gene): 
  data_all <- load_count_data_aggregated_across_celltypes(paths, code_path, opt_aggregation, add_str_pagoda, FALSE, NaN)
    
  #merge different metadata infos
  #impute missing values
  #define columns as factors
  data_all$metadata <- get_metadata(data_all,paths,filename_data_sample_info,filename_data_ancestry, filename_QC_data, opt_load_metadata, opt_scale_and_center_metadata, opt_ancestry_binary)
  paths$results_RUV <-get_results_folder(paths$results_RUV,"","data_split_halves/",NaN,FALSE)
}

#specify RUV factor sequence
if (n_RUV_factors<=10){
  RUV_factors_sequence <- seq(0,n_RUV_factors)
  print(RUV_factors_sequence)
} else{
  RUV_factors_sequence <- c(seq(0,10),seq(15,n_RUV_factors,by=5))
}

#specify design:
model_str_full_ini <- "~ Group"
model_str_red_ini <- "~ 1"

design_formula_full <- formula(eval(parse(text=model_str_full_ini)))
design_formula_red <- formula(eval(parse(text=model_str_red_ini)))

if (opt_run_RUV == TRUE & opt_run_split_halves==TRUE){
  
  ## 1) calculate RUVs and save them
  ## 1a) initial DESeq2 run on data aggregated across all cell types, no RUVs in models (required to calculate RUVs later)

  #remove genes with low counts and store data in DDS objects:
  DDS_all_ini<- get_filtered_data(data_all$CC,data_all$metadata,design_formula_full,paths$xy_chr_list,p_samples_larger_10_counts)
  #determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
  sfe_type <- get_sfe_type(DDS_all_ini)

  flag_DEG_ini <- tryCatch(
    {
      #initially run DSEq2:
      DDS_all_ini <- DESeq(object=DDS_all_ini,parallel=opt_DESeq2_parallel,sfType = sfe_type, test="LRT", reduced = design_formula_red, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
      flag_DEG_ini=0
    },error=function(cond) {
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      flag_DEG_ini <- 303
      return(flag_DEG_ini)
    },
    warning=function(cond) {
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      flag_DEG_ini <- 101
      return(flag_DEG_ini)
    },
    finally={
      message("Initial DEG run is done! \n")
    }
  )
  
  if (flag_DEG_ini<=101){
    ## 1.b) calculate RUVs for n_RUV_factors factors based on filtered data set (including all samples) and save results (to be loaded later)
    res_DDS_ini <- results(DDS_all_ini,alpha=alpha_val,contrast=c("Group","SCZ","CTRL"),independentFiltering = FALSE)
    
    #get data for RUV analysis, only keep genes with p-val >0.5 (clearly non-significant in SCZ/CTRL Group effect testing):
    D <- get_RUV_data(DDS_all_ini,res_DDS_ini,opt_filter_high_CT_specificity_genes,paths$CT_spec_file,filename_CT_spec,FALSE,alpha_val,specificity_TH)
  
    for (n_RUV_factors in RUV_factors_sequence){
      if (n_RUV_factors>0){
        set <- RUVg(D$set_ini, D$empirical, k=n_RUV_factors)
        #save matrix with correlation values after each k_RUV_i result is added:
        saveRDS(set, file=paste0(paths$results_RUV,"set_",as.character(n_RUV_factors),"_RUVS.rds"))
      }
    }
  }
}
if (opt_run_split_halves==TRUE){
  ## 2a) load ct specific data, update paths:
  data_folder <- "CT_clustered_aggregated_data/"
  paths <- get_paths(main_path,main_project_path,data_folder)
  setwd(paths$data)
  countData_filenames <- Sys.glob(paste0("Agg_counts_",add_str_pagoda,"_TH_and_D_adj_filtered",'_*_',opt_aggregation,'.csv'))
  colData_filenames <- Sys.glob(paste0("G_info_",add_str_pagoda,"_TH_and_D_adj_filtered",'_*_',opt_aggregation,'.csv'))
  
  setwd(paths$data)
  data_ct_i <- load_count_data(countData_filenames,colData_filenames,ct_id,opt_aggregation,paths)
  data_ct_i$metadata <- get_metadata(data_ct_i, paths,filename_data_sample_info,filename_data_ancestry, filename_QC_data, opt_load_metadata, opt_scale_and_center_metadata,opt_ancestry_binary)
  
  #create and get results folder:
  celltype <- get_cell_Type_name_from_filename(countData_filenames,ct_id,opt_aggregation)
  paths$results_RUV <-get_results_folder(paths$results_RUV,"","data_split_halves/",NaN,FALSE)
  paths$results_CT_i <-paste0(paths$results_RUV,celltype)
  dir.create(paths$results_CT_i,recursive=TRUE, showWarnings = FALSE)
  
  ## 2b) initial DESeq2 run on data aggregated per cell type, no RUVs in models 
  #remove genes with low counts and located on x and y chromosomes and store data in DDS object:
  DDS_ct_i <- get_filtered_data(data_ct_i$CC,data_ct_i$metadata,design_formula_full,paths$xy_chr_list,p_samples_larger_10_counts)
  #determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
  sfe_type <- get_sfe_type(DDS_ct_i)
  
  #run DSEq2:
  DDS_ct_i <- DESeq(object=DDS_ct_i,
               parallel=opt_DESeq2_parallel,
               sfType = sfe_type, 
               test="LRT", 
               reduced = design_formula_red, 
               useT=TRUE,
               minmu = 1e-6, 
               minReplicatesForReplace = Inf)
  
  #2. for a range of number of RUV factors: run RUV analysis, run deseq2 on random data halves and calculate log2FC
  k_id <- 1
  M_cor_spearman <- array(NaN,c(n_data_splits,length(RUV_factors_sequence)))
  M_cor_pearson <- array(NaN,c(n_data_splits,length(RUV_factors_sequence)))
  for (n_RUV_factors in RUV_factors_sequence){
    if (opt_parallel_for_loop==TRUE){
      cl <- makeCluster(n_workers)
      registerDoParallel(cl)
      M_cor <- foreach(split_i = 1:n_data_splits, .combine='rbind', .errorhandling = 'pass', .packages=c("DESeq2","RUVSeq","doParallel","foreach")) %dopar% {
        X <- get_correlation_LFC_data_halves(DDS_ct_i,
                                             n_RUV_factors,
                                             model_str_full_ini, 
                                             model_str_red_ini, 
                                             opt_DESeq2_parallel,
                                             min_p_donors_groups)
        return(X)
      }
      stopCluster(cl)
      M_cor_spearman[,k_id] <- M_cor[,1]
      M_cor_pearson[,k_id] <- M_cor[,2]
    } else{
      for (split_i in seq(1,n_data_splits)){
        #initialize M_cor:
        M_cor <- array(NaN,c(n_data_splits,2))
  
        colnames(M_cor) <- RUV_factors_sequence
        M_cor[split_i,] <- get_correlation_LFC_data_halves(DDS_ct_i,
                                                            n_RUV_factors,
                                                            model_str_full_ini, 
                                                            model_str_red_ini, 
                                                            opt_DESeq2_parallel,
                                                            min_p_donors_groups)
      }
      M_cor_spearman[,k_id] <- M_cor[,1]
      M_cor_pearson[,k_id] <- M_cor[,2]
    }
  }
  colnames(M_cor_spearman) <- RUV_factors_sequence
  colnames(M_cor_pearson) <- RUV_factors_sequence
  #errors are blocked with parallel option and value in M_cor is set to nan
  #save matrix with correlation values after each k_RUV_i result is added:
  write.table(M_cor[,,1], file=paste0(paths$results_CT_i,"M_cor_spearman.Rdata"))
  write.table(M_cor[,,2], file=paste0(paths$results_CT_i,"M_cor_pearson.Rdata"))
  visualize_RUV_grid_search(M_cor,paths$results_CT_i)

} else {#load results & visualize
  for(cor_method in c("pearson","spearman")){
    M_cor<-read.table(paste0(paths$results_CT_i,"M_cor_",cor_method,".Rdata"))
    visualize_RUV_grid_search(data.matrix(M_cor),paths$results_CT_i)
  }
}



