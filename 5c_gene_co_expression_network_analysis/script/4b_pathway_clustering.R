#about: cluster module GO terms with package rrvgo
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

setwd("../")
main_path = getwd()

setwd(paste0(main_path,"5b_differential_gene_expression_analysis/script/"))
source("pathway_clustering_utils.R")

#get the settings:
opt_input_data = "all_GRN_modules_gprofiler_result"#"GRN_modules_gprofiler_result"
settings <- get_settings_for_input_data(opt_input_data)
p_val_cutoff <- 0.05 # for pathway analysis

run_pathway_clustering_with_rrvgo(opt_input_data, settings, main_path,p_val_cutoff)
