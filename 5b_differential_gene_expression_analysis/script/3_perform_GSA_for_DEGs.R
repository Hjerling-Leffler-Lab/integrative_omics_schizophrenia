#perform GSA analysis based on PFS' updated (2022-07-06) scripts for genes identified with deseq2
#author: Lisa Bast
#date: 15.03.2022
#version: 0.1.5

#load packages:
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)

library(tidyverse)
library(rstatix)

#define paths and load functions:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()

path_DAP_results <- paste0(main_path,"data/proteomics/")
proteomics_hgnc_file <- paste0(main_path,"data/proteomics/hgnc_2022-04-24_small.tsv")
proteomics_filename <- "diagnosis_vs_layers.xlsx" 

set.seed(17274)
graphics.off()

##settings:

version_nr <- "v3"#v2" 

opt_data <- c("transcriptomics") #,"longread",'proteomics') 

alpha_DEGs_DAPs_range <- c(0.3)#this is the best resolution #c(0.05,0.1,0.3,0.5) # for selecting DEGS
alpha_pathways <- 0.05
number_genes_in_gene_set_range <- c(3)#c(3,5,10) 

#for background dataset filter from gene matrix genes that have at least min_meanTPM_BrainCortex expressed in BrainCortex:
min_meanTPM_BrainCortex <- 1 #only used if opt_background_data=="cortex"

options_filter_for_protein_coding_genes <- c(FALSE,TRUE) #

p_val_filter_str <- "P.fdr.group"# "P.bonf.all","P.fdr.all","P.bonf.group","P.fdr.group"

path_background <- paste0(main_path,"/data/") #contains all ensids

for (opt_data_i in opt_data){
  if (opt_data_i == 'transcriptomics'){
    opt_background_data_sets <- c("cortex")#,"SnRNAseqDet") 
    add_folder_str <- "DEGs/"
  } else if (opt_data_i == "proteomics"){
    opt_create_background_from_proteomics <- TRUE #FALSE #
    proteomics_data <- load_proteomics_data(path_DAP_results,proteomics_filename, proteomics_hgnc_file)
    #add ensgid to proteomics_data
    #create background out of all genes mapping to proteins tested:
    if (opt_create_background_from_proteomics){
      create_background_from_proteomics(proteomics_data,path_background)
    }
    opt_background_data_sets <- c("cortex")#,"genesMappingToProteinsDetectedInProteomics") 
    add_folder_str <- "DAPs/"
  } else if (opt_data_i == "longread"){
    opt_background_data_sets <- c("cortex")
    path_lr_results <- paste0(main_path,"output/DIUGs/")
    filename <- "sr_mm_PC12_filtE_gene_FDR.tsv"
    lr_data <- load_longread_data(path_lr_results,filename)
    add_folder_str <- "DIUGs/"
  }
}

for (opt_filter_for_protein_coding_genes in options_filter_for_protein_coding_genes){
  for (opt_background_data in opt_background_data_sets){
    if (opt_filter_for_protein_coding_genes==TRUE){
      add_str_background = paste0(paste0('/',opt_background_data),"_pc")
      genetype='protein_coding'
    }else{
      add_str_background = paste0('/',opt_background_data)
      genetype='all'
    }
    add_str_velo <- "_cellranger"
    n_cluster_vec <- c(16)#,37,3)
    
    n_RUV_vars_to_inlcude_range <- c(6) #c(0,4,6,8) 
    
    #paths to data files:
    path_geneset <- paste0(main_path,"/data/","gsa-files-",version_nr,"/") #contains info about go terms, paths and ensgid (GO term database)
    
    for (alpha_val in alpha_DEGs_DAPs_range){
      #filename definition:
      if (opt_data_i == 'transcriptomics'){
        filename <- "df_p_vals_adj_sum_per_celltype.csv"
      } 
      
      #load genesets and geneMatrix file
      #unused at the moment: 
      #gene_sets_meta <- LOAD("genesets.metadata.tsv",path_geneset)
      if (version_nr=='v3'){
        gene_sets_str <- "gene_setsA, gene_setsB,"
        gene_setsA <- LOAD("genesetsA.tsv",path_geneset)
        gene_setsB <- LOAD("genesetsB.tsv",path_geneset)
        #make sure each GO category is a separate group:
        gene_setsA$group[gene_setsA$group=="gene-ontology"] <- paste(gene_setsA$group[gene_setsA$group=="gene-ontology"],gene_setsA$subgroup[gene_setsA$group=="gene-ontology"],sep='_')
        gene_setsB$group[gene_setsB$group=="gene-ontology"] <- paste(gene_setsB$group[gene_setsB$group=="gene-ontology"],gene_setsB$subgroup[gene_setsB$group=="gene-ontology"],sep='_')
      } else{
        gene_sets_str <- "gene_sets,"
        gene_sets <- LOAD("genesets.tsv",path_geneset)
        gene_sets$group[gene_sets$group=="gene-ontology"] <- paste(gene_sets$group[gene_sets$group=="gene-ontology"],gene_sets$subgroup[gene_sets$group=="gene-ontology"],sep='_')
      }
      
      background <- get_background(path_background,min_meanTPM_BrainCortex,genetype,opt_background_data,opt_filter_for_protein_coding_genes,version_nr)

      if (opt_data_i == 'transcriptomics'){
        str_true <- "True"
        str_false <- "False"
        #perform GSA for DESeq2 DEGs:
        for (n_cluster in n_cluster_vec){
          for (min_number_genes_in_gene_set in number_genes_in_gene_set_range){
            for (n_RUV_vars_to_inlcude in n_RUV_vars_to_inlcude_range){
              if (n_RUV_vars_to_inlcude<=8){ # otherwise overfitting
                bool_first <- data.frame (up = c(TRUE),down = c(TRUE))
                #define design str:
                if (n_RUV_vars_to_inlcude >0){
                  add_str_design <- paste0('design_with_',n_RUV_vars_to_inlcude,"_RUVs")
                  add_str_design_short <- paste0(n_RUV_vars_to_inlcude,"_RUVs")
                } else{
                  add_str_design <- "design_without_covariates"
                  add_str_design_short <- "no_covariates"
                }
                #DEG results path and GSA results_path definition:
                if (n_cluster==3){
                  path_DEG_results <- paste0(main_path,"/output/DEGs/3_classes/",add_str_design,"/")
                  path_GSA_results <- paste0(main_path,"/output/GSA_analysis/",add_folder_str,version_nr,"/3_classes/alpha",sub("\\.", "", as.character(alpha_val*10)),"/min_overlap_",as.character(min_number_genes_in_gene_set),"/",add_str_design_short,add_str_background,"/")
                } else{
                  path_DEG_results <- paste0(main_path,"/output/DEGs/",n_cluster,"_CTs/",add_str_design,"/")
                  path_GSA_results <- paste0(main_path,"/output/GSA_analysis/",add_folder_str,version_nr,'/',n_cluster,"_CTs/alpha",sub("\\.", "", as.character(alpha_val*10)),"/min_overlap_",as.character(min_number_genes_in_gene_set),"/",add_str_design_short,add_str_background,"/")
                }
                
                print(path_GSA_results)
                if (dir.exists(path_GSA_results)==FALSE){
                  dir.create(path_GSA_results,recursive=TRUE) 
                }
                
                #load DEG results and extract ensgids of significant genes
                df_DEG = read.csv(paste0(path_DEG_results,filename))
                print(paste0("number of DEGs: ",as.character(length(unique(df_DEG$Gene_short)))))
                #rename columns, replace dot with underscore
                for (column_i in colnames(df_DEG)){
                  colnames(df_DEG)[colnames(df_DEG)==column_i] <- gsub("\\.",'_',column_i)
                }
                
                #bool columns have no value if nan --> problematic?
                pattern = c("^Exc","^Inh","^Non","^Astro","^Endo","^Micro","^Oligo")
                CT_cluster_list <- c()
                for (p in pattern){
                  CT_cluster_list <- c(CT_cluster_list,colnames(df_DEG) %>% str_subset(pattern=p))
                }
        
                #for each cell type cluster perform GSA:
                for (CT_str in CT_cluster_list){
                  #add column indicating significance of up-regulated genes:
                  
                  df_DEG_up <- eval(parse(text=paste0("df_DEG %>% mutate(",CT_str,"_up_significant"," = ifelse((",CT_str," <= alpha_val) & (bool_",CT_str,"_upregulated %in% str_true), TRUE, FALSE)",")")))
                  #add column indicating significance of down-regulated genes:
                  df_DEG_down <- eval(parse(text=paste0("df_DEG %>% mutate(",CT_str,"_down_significant"," = ifelse((",CT_str," <= alpha_val) & (bool_",CT_str,"_upregulated %in% str_false), TRUE, FALSE)",")")))
                  #add column indicating significance of de-regulated (up- or down-reg.) genes:
                  df_DEG <- eval(parse(text=paste0("df_DEG %>% mutate(",CT_str,"_significant"," = ifelse((",CT_str," <= alpha_val), TRUE, FALSE)",")")))
                  
                  #  df_DEG <- df_DEG %>% mutate(Inhibitory_significant = if_else(Inhibitory <= alpha_val, TRUE, FALSE))
                  df_DEG_CT_up <- df_DEG_up %>% select("ensgid",paste0(CT_str,"_up_significant"))
                  df_DEG_CT_down <- df_DEG_down %>% select("ensgid",paste0(CT_str,"_down_significant"))
                  df_DEG_CT <- df_DEG %>% select("ensgid",paste0(CT_str,"_significant"))
                  
                  #remove all genes from DEG list that are not part of the background:
                  #this indirectly filters for protein coding genes if opt_filter_for_protein_coding_genes is set to TRUE
                  df_DEG_CT_up <- subset(df_DEG_CT_up, ensgid %in% background$ensgid)
                  df_DEG_CT_down <- subset(df_DEG_CT_down, ensgid %in% background$ensgid)
                  df_DEG_CT <- subset(df_DEG_CT, ensgid %in% background$ensgid)
                  
                  if (dim(df_DEG_CT_up)[1]>1){
                    #add df_DEG_CT lists to background matrix as column:
                    background <- full_join(background, df_DEG_CT_up,by="ensgid")
                  }
                  if (dim(df_DEG_CT_down)[1]>1){
                    background <- full_join(background, df_DEG_CT_down,by="ensgid")
                  }
                  if (dim(df_DEG_CT)[1]>1){
                    background <- full_join(background, df_DEG_CT,by="ensgid")
                  }
                  
                  if (eval(parse(text=paste0("sum(background$",CT_str,"_up_significant,na.rm=TRUE)>0")))){
                    # 1) GSA of transcription factor genes
                    #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
                    result_GSA_pathways_up <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",CT_str,"_up_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
                    # select pathways
                    # there are a LOT of pathways, pick those of interest and not redundant 
                    result_GSA_sign_pathways_up <- select_pathways(result_GSA_pathways_up,alpha_pathways,version_nr,p_val_filter_str)
                    print(CT_str)
                    print(' upregulated: ')
                    print(paste0(dim(result_GSA_sign_pathways_up)[1]," out of ",dim(result_GSA_pathways_up)[1]))
                    
                    if (dim(result_GSA_pathways_up)[1]>0){
                      save_pathways_as_csv_file(result_GSA_pathways_up,result_GSA_sign_pathways_up, path_GSA_results,CT_str,'up')
                    }
                  }
                  
                  if (eval(parse(text=paste0(paste0("sum(background$",CT_str),"_down_significant,na.rm=TRUE)>0")))){
                    # 1) GSA of transcription factor genes
                    #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
                    result_GSA_pathways_down <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",CT_str,"_down_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
                    # select pathways
                    # there are a LOT of pathways, pick those of interest and not redundant 
                    result_GSA_sign_pathways_down <- select_pathways(result_GSA_pathways_down,alpha_pathways,version_nr,p_val_filter_str)
                    print(CT_str)
                    print(' downregulated: ')
                    print(paste0(dim(result_GSA_sign_pathways_down)[1]," out of ",dim(result_GSA_pathways_down)[1]))
                    
                    if (dim(result_GSA_pathways_down)[1]>0){
                      save_pathways_as_csv_file(result_GSA_pathways_down,result_GSA_sign_pathways_down, path_GSA_results,CT_str,'down')
                    }
                  }
                  
                  if (eval(parse(text=paste0(paste0("sum(background$",CT_str),"_significant,na.rm=TRUE)>0")))){
                    # 1) GSA of transcription factor genes
                    #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
                    result_GSA_pathways <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",CT_str,"_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
                    # select pathways
                    # there are a LOT of pathways, pick those of interest and not redundant 
                    result_GSA_sign_pathways <- select_pathways(result_GSA_pathways,alpha_pathways,version_nr,p_val_filter_str)
                    print(CT_str)
                    print(' upregulated: ')
                    print(paste0(dim(result_GSA_sign_pathways)[1]," out of ",dim(result_GSA_pathways)[1]))
                    
                    if (dim(result_GSA_pathways)[1]>0){
                      save_pathways_as_csv_file(result_GSA_pathways,result_GSA_sign_pathways, path_GSA_results,CT_str,'')
                    }
                  }
                }
              }
            }
          }
        }
      } 
      for (min_number_genes_in_gene_set in number_genes_in_gene_set_range){
        if (opt_data_i == "proteomics"){
          path_GSA_results <- paste0(main_path,"output/GSA_analysis/",add_folder_str,version_nr,"/alpha",sub("\\.", "", as.character(alpha_val*10)),"/min_overlap_",min_number_genes_in_gene_set,"/",add_str_background,"/")
          print(path_GSA_results)
          if (dir.exists(path_GSA_results)==FALSE){
            dir.create(path_GSA_results,recursive=TRUE) 
          }
          
          for (m in c('both','down','up')){
            if (m=='down'){
              proteomics_data_sel <- filter(proteomics_data, log2fc<0)
            } else if (m=='up'){
              proteomics_data_sel <- filter(proteomics_data, log2fc>0)
            } else{
              proteomics_data_sel <- proteomics_data
            }

            #df_DAP = proteomics_data_sel %>% select(protein_id, symbol = gene_symbols_or_id, hgnc_id, ensgid, peptides_used_for_dea, log2fc, qvalue)
            df_DAP = proteomics_data_sel %>% select(protein_id, symbol = hgnc_symbol, hgnc_id, ensgid, peptides_used_for_dea, log2fc, qvalue)
            
            #add column indicating significance of up-/ down-regulated genes:
            df_DAP <- eval(parse(text=paste0("df_DAP %>% mutate(",m,"_significant = ifelse((qvalue <=",as.character(alpha_val),"),TRUE,FALSE))")))
            
            #  df_DEG <- df_DEG %>% mutate(Inhibitory_significant = if_else(Inhibitory <= alpha_val, TRUE, FALSE))
            df_DAP_sel <- df_DAP %>% select("ensgid",paste0(m,"_significant"))
            
            #remove all genes from DEG list that are not part of the background:
            #this indirectly filters for protein coding genes if opt_filter_for_protein_coding_genes is set to TRUE
            df_DAP_sel <- subset(df_DAP, ensgid %in% background$ensgid)
            
            #add df_DEG_CT lists to background matrix as column:
            background <- full_join(background, df_DAP_sel,by="ensgid")
            
            if (eval(parse(text=paste0(paste0("sum(background$",m),"_significant,na.rm=TRUE)>0")))){
              # 1) GSA of transcription factor genes
              #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
              result_GSA_pathways <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",m,"_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
              # select pathways
              # there are a LOT of pathways, pick those of interest and not redundant 
              result_GSA_sign_pathways <- select_pathways(result_GSA_pathways,alpha_pathways,version_nr,p_val_filter_str)
              print(paste0(m,'-regulated: '))
              print(paste0(paste0(dim(result_GSA_sign_pathways)[1]," out of "),dim(result_GSA_pathways)[1]))
              
              #save results in excel workbooks
              if (dim(result_GSA_pathways)[1]>0){
                save_pathways_as_csv_file(result_GSA_pathways,result_GSA_sign_pathways, path_GSA_results,'proteomics',m)
              }
            }
          }
        } else if (opt_data_i =="longread"){
            m="both"
            
            path_GSA_results <- paste0(main_path,"output/GSA_analysis/",add_folder_str,version_nr,"/alpha",sub("\\.", "", as.character(alpha_val*10)),"/min_overlap_",min_number_genes_in_gene_set,add_str_background,"/")
            print(path_GSA_results)
            if (dir.exists(path_GSA_results)==FALSE){
              dir.create(path_GSA_results,recursive=TRUE) 
            }
            
            #df_DAP = proteomics_data_sel %>% select(protein_id, symbol = gene_symbols_or_id, hgnc_id, ensgid, peptides_used_for_dea, log2fc, qvalue)
            df_lr = lr_data %>% select(associated_gene,FDR)
            
            #add column indicating significance of up-/ down-regulated genes:
            df_lr <- eval(parse(text=paste0("df_lr %>% mutate(",m,"_significant = ifelse((FDR <=",as.character(alpha_val),"),TRUE,FALSE))")))
            
            #  df_DEG <- df_DEG %>% mutate(Inhibitory_significant = if_else(Inhibitory <= alpha_val, TRUE, FALSE))
            df_lr_sel <- df_lr%>% select("associated_gene",paste0(m,"_significant"))
            
            #remove all genes from DEG list that are not part of the background:
            #this indirectly filters for protein coding genes if opt_filter_for_protein_coding_genes is set to TRUE
            
            df_lr_sel <- subset(df_lr, associated_gene %in% background$gene_name)
            
            #rename column
            names(df_lr_sel)[names(df_lr_sel) == "associated_gene"] <- "gene_name"
            
            #add df_DEG_CT lists to background matrix as column:
            background <- full_join(background, df_lr_sel,by="gene_name")
            
            if (eval(parse(text=paste0(paste0("sum(background$",m),"_significant,na.rm=TRUE)>0")))){
              # 1) GSA of transcription factor genes
              #b$Excitatory is a logical column, TRUE=query genes for gene set analysis that indicates which genes we want to test for
              result_GSA_pathways <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",m,"_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
              #result_GSA_pathways <- eval(parse(text=paste0(" GSA_",version_nr,"(background, ", "\"",m,"_significant\"",", ",gene_sets_str,as.character(min_number_genes_in_gene_set),")"))) 
              
              # select pathways
              # there are a LOT of pathways, pick those of interest and not redundant 
              result_GSA_sign_pathways <- select_pathways(result_GSA_pathways,alpha_pathways,version_nr,p_val_filter_str)
              print(paste0('differential isoform usage pathways: '))
              print(paste0(dim(result_GSA_sign_pathways)[1]," out of ",dim(result_GSA_pathways)[1]))
              
              #save results in excel workbooks
              if (dim(result_GSA_pathways)[1]>0){
                save_pathways_as_csv_file(result_GSA_pathways,result_GSA_sign_pathways, path_GSA_results,'longread',m)
              }
            }
          }
      }
    }
  }
}



