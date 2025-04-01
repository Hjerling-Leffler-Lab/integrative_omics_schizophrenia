#about: load result from GRN inference
#       perform GSA of significant gene sets 
#author: Lisa Bast
#date: 08.03.2022
#version: 0.0.2

library(scITD)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(sjmisc)

library(tidyverse)
library(rstatix)

version_nr <- "v3"

set.seed(17274)
graphics.off()

path_code = getwd()
setwd("../")
main_path_project  = getwd()
setwd("../")
main_path  = getwd()

path_genelist_to_test <- paste0(main_path_project,"/output/cluster_name_15CTs/scz/cor_bicor/")
results_path <- paste0(path_genelist_to_test,"/GSA_results/")
#make sure results path exists
if (file.exists(results_path)==FALSE) {
  dir.create(results_path, recursive=TRUE)
}

#load helper functions
source(paste0(path_code,'utils_hdWGCNA.R'))
source(paste0(main_path,"/5b_differential_gene_expression_analysis/script/",'utils.R'))

#define settings:
n_cores <- 12

#GSA settings
path_background <- paste0(main_path,"data/Gene_lists/") #contains all ensids
opt_background_data <- "BrainCortexExpressed"
number_genes_in_gene_set_range <- c(3)#c(3,5,10) 
alpha_pathways <- 0.05
p_val_filter_str <- "P.fdr.group"# "P.bonf.all","P.fdr.all","P.bonf.group","P.fdr.group"

#for background dataset filter from gene matrix genes that have at least min_meanTPM_BrainCortex expressed in BrainCortex:
min_meanTPM_BrainCortex <- 1 #only used if opt_background_data=="BrainCortexExpressed"
options_filter_for_protein_coding_genes <- c(FALSE) #c(TRUE,FALSE)

for (opt_filter_for_protein_coding_genes in options_filter_for_protein_coding_genes){
  if (opt_filter_for_protein_coding_genes==TRUE){
    add_str_background = paste0(paste0('/',opt_background_data),"_protein_coding_only")
    genetype='protein_coding'
  } else{
    add_str_background = paste0('/',opt_background_data)
    genetype='all'
  }
  
  #paths to data files:
  path_geneset <- paste0(path_code,"gsa-files-",version_nr,"/") #contains info about go terms, paths and ensgid (GO term database)
  
  #load genesets and geneMatrix file
  #unused at the moment: 
  #gene_sets_meta <- LOAD("genesets.metadata.tsv",path_geneset)
  if (version_nr=='v3'){
    gene_sets_str <- "gene_setsA, gene_setsB,"
    gene_setsA <- LOAD("genesetsA.tsv",path_geneset)
    gene_setsB <- LOAD("genesetsB.tsv",path_geneset)
    #make sure each GO category is a separate group:
    #gene_setsA$group[gene_setsA$group=="gene-ontology"] <- paste(gene_setsA$group[gene_setsA$group=="gene-ontology"],gene_setsA$subgroup[gene_setsA$group=="gene-ontology"],sep='_')
    gene_setsB$group[gene_setsB$group=="gene-ontology"] <- paste(gene_setsB$group[gene_setsB$group=="gene-ontology"],gene_setsB$subgroup[gene_setsB$group=="gene-ontology"],sep='_')
  } else{
    gene_sets_str <- "gene_sets,"
    gene_sets <- LOAD("genesets.tsv",path_geneset)
    gene_sets$group[gene_sets$group=="gene-ontology"] <- paste(gene_sets$group[gene_sets$group=="gene-ontology"],gene_sets$subgroup[gene_sets$group=="gene-ontology"],sep='_')
  }
  background <- get_background(path_background,min_meanTPM_BrainCortex,genetype,opt_background_data,opt_filter_for_protein_coding_genes,version_nr)
  
  str_true <- "True"
  str_false <- "False"
  
  for (min_number_genes_in_gene_set in number_genes_in_gene_set_range){

    gene_list_to_test <- get_gene_lists_to_test(path_genelist_to_test)
    #remove all genes from DEG list that are not part of the background:
    #this indirectly filters for protein coding genes if opt_filter_for_protein_coding_genes is set to TRUE
    gene_list_to_test <- subset(gene_list_to_test, ensgid %in% background$ensgid)
    if (dim(gene_list_to_test)[1]>1){
      #add gene lists to background matrix as column:
      background_tmp <- full_join(background, gene_list_to_test,by="ensgid")
      intersect_strings <- colnames(gene_list_to_test)
      intersect_strings <- intersect_strings[intersect_strings!="ensgid"]
      for (intersect_str in intersect_strings){
        if (eval(parse(text=paste0("sum(background_tmp$",intersect_str,",na.rm=TRUE)>0")))){
          # 1) GSA of transcription factor genes
          #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
          
          result_GSA_pathways <- eval(parse(text=paste0(" GSA_",version_nr,"(background_tmp, ", "\"",intersect_str,"\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
          # select pathways
          # there are a LOT of pathways, pick those of interest and not redundant 
          result_GSA_sign_pathways <- select_pathways(result_GSA_pathways,alpha_pathways,version_nr,p_val_filter_str)
          print(paste0(intersect_str,': '))
          print(paste0(paste0(dim(result_GSA_sign_pathways)[1]," out of "),dim(result_GSA_pathways)[1]))
          #browser()
          if (dim(result_GSA_pathways)[1]>0){
            save_pathways_as_csv_file(result_GSA_pathways,result_GSA_sign_pathways, results_path, intersect_str,"")
          }
        }
      }
      rm(gene_list_to_test)
      rm(background_tmp)
      rm(result_GSA_pathways)
      rm(result_GSA_sign_pathways)
    }
  }
}

