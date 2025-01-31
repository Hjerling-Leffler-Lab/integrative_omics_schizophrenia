#author: Lisa Bast
#date: 13.06.2023
#about: export ensgids for (inner/ outer) mitochondria, neuron differentiation and immune response GO terms
#version: 0.0.1

library(data.table)
library(tidyverse)
library(readxl)
library(clusterProfiler)

#define paths and load functions:
code_path = getwd()
source("utils.R")

setwd("../")
main_path_project = getwd()

path.out <- paste0(main_path_project,"/Gene_lists_GO_DB/")
GO_CC_filename <- paste0(path.out,"hsapiens.GO_CC.ENSG.gmt")
GO_BP_filename <- paste0(path.out,"hsapiens.GO_BP.ENSG.gmt")

df_CC <- read.gmt(GO_CC_filename) 
df_BP <- read.gmt(GO_BP_filename) 

search_terms_CC <- c("mitochondrial outer membrane","mitochondrial inner membrane")
search_term_IDs_CC <- c("GO:0005741","GO:0005743")

search_terms_BP <- c("neuron differentiation","immune response")#
search_term_IDs_BP <- c("GO:0030182","GO:0006955")

save_GO_ensgids(df_CC,search_terms_CC,search_term_IDs_CC,"CC")
save_GO_ensgids(df_BP,search_terms_BP,search_term_IDs_BP,"BP")
