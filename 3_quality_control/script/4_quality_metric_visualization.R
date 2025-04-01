
#author: Lisa Bast
#date: 12.10.2020
#about visualization of quality metrics for SCz-CTRL snRNAseq data set, Fig. S1 E,F
#version: 0.0.2

library(GGally)
library(corrgram)
#library(ellipse)
library(corrplot)
library(factoextra) #pca visualization

rm(list=ls())
opt_without_outlier=TRUE

#define paths
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()

path_code <- paste0(main_project_path,"/3_quality_control/script/")
path_results <- paste0(main_project_path,"/3_quality_control/output/")


#read metrics per donor sample:
metrics <- read.csv(file = paste0(path_results,'metrics_summary_tidy.csv'))
          
selected_cols <- which(names(metrics)%in%c("num_umi",
                                           "mean_reads_per_umi",
                                           "num_reads",
                                           "p_valid_barcodes",
                                           "p_sequencing_saturation",
                                           "p_genome_not_gene",
                                           "p_mapped_reads",
                                           "p_unmapped_reads",
                                           "p_cell_del_filt",
                                           "mean_counts_per_barcode",
                                           "std_counts_per_barcode",
                                           "median_gpc" ))

if (opt_without_outlier==TRUE) {
  #selected_rows <- metrics$donor_ID_python!='S66' & metrics$donor_ID_python!='S78'
  metrics <- subset(subset(metrics, donor_ID_python!='S66'), donor_ID_python!='S78')
  fig_str <- "_without outlier.png"
}else{
  fig_str <-".png"
}

#graphics:
plot_correlation_of_metrics(metrics, selected_cols,path_results,fig_str,"spearman",FALSE)
plot_correlation_of_metrics(metrics, selected_cols,path_results,fig_str,"spearman",TRUE)
plot_correlation_of_metrics(metrics, selected_cols,path_results,fig_str,"pearson",FALSE)
plot_pearson_vs_spearman_correlation(metrics,selected_cols,path_results,fig_str)
plot_PCA_of_samples(metrics,selected_cols,path_results,TRUE,TRUE,"Group")


