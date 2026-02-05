#about: cluster GO terms with package rrvgo
#       Fig. 3A,C, S8A-C
#author:Lisa Bast
#date: 22.09.2022
#version: 0.0.3

#load libraries
library(rrvgo)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
#library(GOSemSim)

#define paths and load functions:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()

#get the settings:
opt_input_data = "DEG_GSA_result_all_CTs"#"DEG_GSA_result"#"all_GRN_modules_gprofiler_result"#"GRN_modules_gprofiler_result"##"GRN_modules_GSA_result"#"all_GRN_modules_GSA_result" #  ##"longread_results"###"DEG_GSA_and_GRN_modules_gprofiler_result"####"Patient_clustering_genes_GSA_result"##"Patient_clustering_genes_GSEA_result"#
settings <- get_settings(opt_input_data)
p_val_cutoff <- 0.05 # for pathway analysis

run_pathway_clustering_with_rrvgo(opt_input_data, settings, main_path,p_val_cutoff)