#author: Lisa Bast
#date: 13.06.2023
#about: export ensgids for (inner/ outer) mitochondria GO terms
#version: 0.0.1

library(data.table)
library(tidyverse)
library(readxl)
library(clusterProfiler)

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

path_code = getwd()
setwd("../")
main_path_project  = getwd()

path.out <- paste0(main_path_project,"/output/Gene_lists/GO_DB/")
GO_CC_filename <- paste0(path.out,"hsapiens.GO_CC.ENSG.gmt")

search_terms_CC <- c("mitochondrion")
search_term_IDs_CC <- c("GO:0005739")

df_CC <- read.gmt(GO_CC_filename) 

save_GO_ensgids(df_CC,search_terms_CC,search_term_IDs_CC,"CC")

