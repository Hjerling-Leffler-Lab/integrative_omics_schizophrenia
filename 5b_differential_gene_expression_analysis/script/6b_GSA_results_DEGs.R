#---
#author: Shuyang Yao, modified by Lisa Bast
#output: enrichment bubble plot, Fig. 3A
#date: "2025-05-09"
#---


#- read in DEG GSA result for alpha value
#- plot for BP, CC, MF separately, divide by regulation (up/ down)
#- optionally run rrvgo to structure GO terms

library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

#define paths and load functions:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()

# set paths
DEG.alphas=c(0.05,0.1,0.3)
GSA.FDRthresh=0.05
GO_groups <- c("BP","CC","MF")
opt_load_rrvgo <- TRUE#FALSE#
if (opt_load_rrvgo==FALSE){
  opt_sim_based_on = "scores" #'size' #
  rrvgo_ths <- c(0.7,0.8,0.9)
} else{
  opt_sim_based_on = "scores"
  rrvgo_ths <- c('loaded')
}


for (DEG.alpha in DEG.alphas){
  if (DEG.alpha<0.1){
    WD = paste0(main_path,"/output/GSA_analysis/DEGs/v3/alpha",as.character(gsub("0.0","0",DEG.alpha)),"/min_overlap_3/6_RUVs/cortex/")
  } else{
    WD = paste0(main_path,"/output/GSA_analysis/DEGs/v3/alpha",as.character(gsub("0.","",DEG.alpha)),"/min_overlap_3/6_RUVs/cortex/")
  }
  for (GO_group in GO_groups){
    for (rrvgo_th in rrvgo_ths){
      dat <- load_data(WD,GO_group,GSA.FDRthresh,opt_load_rrvgo,opt_sim_based_on,rrvgo_th)
      # columns of data, only list columns useful for the plot
      colnames(dat)
      #- "celltype_name" (cell type, e.g., Ex_3-4_IT)        
      #- "geneset" (e.g., GO:0099504 synaptic vesicle cycle)
      #- "P.fdr.group" (FDR based on P values, got from our standard GSA script)
      #add number of genes in GO term
      #- "nGenes_GO" (number of genes in the GO term, got from our standard GSA script)
      #- "nGenes_CT" (number of genes in the cell type, got from our standard GSA script)
      #- "nGenes_overlap" (number of genes in the intersetion, got from our standard GSA script)
      #- "parentTerm" (name of parent GO term, got by running rrvgo on GSA results.)
      #- "nGenes_overlap_pct" (derived from the above nGenes_* columns; = nGenes_overlap/(nGenes_CT + nGenes_GO - nGenes_overlap))
      
      # plot: 
      plot_GSA_enrichment(dat,WD,DEG.alpha,GSA.FDRthresh,GO_group,T,opt_sim_based_on,rrvgo_th,"DEGs")
      if (DEG.alpha==0.05){
        plot_GSA_enrichment(dat,WD,DEG.alpha,GSA.FDRthresh,GO_group,F,opt_sim_based_on,rrvgo_th,"DEGs")
      }
    }
  }
}

