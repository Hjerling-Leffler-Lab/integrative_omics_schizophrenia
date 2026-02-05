#about: calculate RUV factors for specific number of RUV factors
#       on whole data
#       prior to DESeq2 run
#author: Lisa Bast
#date: 15.02.22
#version: 0.0.1

library('DESeq2')
library("RUVSeq")#RUV analysis
library("variancePartition")#to analyse each factors contribution to the expression variation of each gene
library('readxl')
library("tidyverse")#for transforming matrix into tibble
library("foreach")#parallel for loops
library("doParallel")
library("stringr")
library("iotools") # requires for read.csv.raw efficiently read csv files
library("corrplot")
library("beeswarm")
library("plotly")
library('mice')
#library("gridExtra")
#library("ggplot2")
#library("ggbeeswarm")

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()

graphics.off()

n_workers <- 6 #20 #35

n_cluster <- 15
opt_pagoda <- TRUE#FALSE#
if (opt_pagoda==TRUE){
  add_str_pagoda <- '_pagoda'
} else{
  add_str_pagoda <- ''
}

opt_load_metadata <- FALSE# #TRUE# otherwise missing values are imputed 
opt_ancestry_binary <- TRUE

#RUV settings:
opt_filter_high_CT_specificity_genes <- TRUE
opt_run_RUV = TRUE#FALSE#
pval_TH <- 0.5
specificity_TH <- 0.5 #0.8

#DESEQ2 settings:
alpha_val <- 0.05

#RUV settings:
n_RUV_factors <- 6 

#which cell type annotation files should be used?:
opt_celltype_groups = 'conos_cluster_based'#'scmap'
opt_aggregation = "sum" #'mean'

if (opt_aggregation == 'sum'){
  opt_filter_genes_with_low_counts <- TRUE
} else{
  opt_filter_genes_with_low_counts <- FALSE
}
if (opt_filter_genes_with_low_counts){
  p_samples_larger_10_counts<-70
}
opt_parallel_for_loop = TRUE # for RUV grid search 
opt_DESeq2_parallel = FALSE # causes some errors I don't understand. better set to FALSE for now
opt_scale_and_center_metadata = TRUE
opt_plot_variance = TRUE
#which metadata files should be used?:
filename_data_sample_info <- "T1_basic_donor_information.xlsx" 
filename_QC_data <- "T2_quality_metrics_per_sample.xlsx"  
filename_data_ancestry <- 'sample_brn_ancestry.tsv'
filename_CT_spec <- 'speMx_Expr_15CTs.Rdata'

#main_subfolder <- "5b_differential_gene_expression_analysis"
data_folder <- "Aggregated_data/"
opt_downsampled_data == FALSE
subfolder <- ""
i <- NaN

paths <- get_paths(main_path,main_project_path,data_folder)


## 0) load data:
# a) expression data (aggregated over all cells from a specific sample for each gene): 
data_all <- load_count_data_aggregated_across_celltypes(paths, code_path, opt_aggregation, add_str_pagoda, opt_downsampled_data, i)


# b) metadata:
#merge different metadata infos
#impute missing values
#define columns as factors
data_all$metadata <- get_metadata(data_all,paths,filename_data_sample_info,filename_data_ancestry, filename_QC_data, opt_load_metadata, opt_scale_and_center_metadata,opt_ancestry_binary)


## initial deseq2 run without any RUVs to determine genes with large p-values (p-val>0.5) and select them for RUV estimation
#a) load data in deseq2 object, remove genes with low counts (at least p% of samples need to have more than 10 counts), remove genes on xy chromosomes and store data in DDS objects:
DDS<- get_filtered_data(data_all$CC,data_all$metadata,"~ Group",paths$xy_chr_list,opt_filter_genes_with_low_counts)

#b) determine how many genes have positive counts (>0) for all samples and decide for sfe setting based on that
sfe_type <- get_sfe_type(DDS)

#c) run  deseq2
flag_DEG_ini <- tryCatch(
  {
    #initially run DSEq2:
    DDS <- DESeq(object=DDS,parallel=opt_DESeq2_parallel,sfType = sfe_type, test="LRT", reduced = ~1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
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
#d) obtain DEseq2 results (for initial run)
if (flag_DEG_ini<=101){
  res_RUV <- results(DDS,alpha=alpha_val,contrast=c("Group","SCZ","CTRL"),independentFiltering = FALSE)
}

#e) get data for RUV analysis, only keep genes with p-val >0.5 (clearly non-significant in SCZ/CTRL Group effect testing):
D <- get_RUV_data(DDS,res_RUV,opt_filter_high_CT_specificity_genes,paths$CT_spec_file,filename_CT_spec,FALSE,pval_TH,specificity_TH)

## perform RUV analysis on whole data set for range of RUV factors and store result
#create and get results subfolder:
design_interactions_str <- paste0('design_with_',as.character(n_RUV_factors),'RUV_covariates_',opt_aggregation)
paths$res_RUV <- get_results_folder(paths$results_RUV,design_interactions_str, data_all$results_subfolder,i,FALSE)
#make sure path exists
dir.create(paths$res_RUV, recursive=TRUE, showWarnings = FALSE)
RUV_factors <- calculate_and_export_RUVs(D,DDS,data_all$colData,n_RUV_factors,paths$res_RUV, opt_plot_variance, opt_downsampled_data ,as.integer(i))

