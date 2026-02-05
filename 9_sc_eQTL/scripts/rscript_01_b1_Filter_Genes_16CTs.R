# File name: 00_01_b1_Filter_Genes_15CTs.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-04-02 
# Aim: Remove genes for eQTL analysis; keep genes with ≥1 count in ≥10% cells in each cell type

#--- script starts ---#

# libraries 
library(tidyverse)
library(loomR)
library(data.table)

# specify paths
filepath <- "" # folder of loom files
filelist <- grep(list.files(filepath, pattern = "loom"), pattern='downsampled', invert=TRUE, value=TRUE)
WD <- "./9_sc_eQTL/data/eQTL/pheno/16CT" # working directory
setwd(WD)

# read in loom files
# input loom files
for(i in 1:length(filelist)){
  file=paste(filepath, filelist[i], sep="")
  mycluster=basename(file)
  mycluster <- gsub("Samples_cellranger_pagoda_TH_and_D_adj_filtered_conos_cluster_based_","",mycluster)
  mycluster <- gsub(".loom","",mycluster)
  
  # link to loom file
  lfile <- connect(file,mode = "r", skip.validate = TRUE)
  # get expression_cell matrix
  exp <- as.data.frame(t(lfile[["matrix"]][,]))
  # number of zeros per row
  exp2 <- rowSums(exp==0)%>% as.data.frame()
  colnames(exp2) <- "count_cell_0exp"
  exp2$N_cell_total <- ncol(exp)
  rm(exp)
  exp2$ENSGID <- lfile[["row_attrs/Accession"]][]
  exp2$cluster <- mycluster
  # combine
  if(i==1){
    mydf <- exp2
    rm(exp2)
    print(i)
  }else{
    mydf <- rbind(mydf, exp2)
    rm(exp2)
    print(i)
  }
}

# per cell type, filter out: genes with count_cell_0exp > 0.9*n_total_cells
# output: genes to keep per cell type

gene2keep <- mydf %>%
  filter(count_cell_0exp<=0.9*N_cell_total)

for(i in 1:length(filelist)){
  file=paste(filepath, filelist[i], sep="")
  mycluster=basename(file)
  mycluster <- gsub("Samples_cellranger_pagoda_TH_and_D_adj_filtered_conos_cluster_based_","",mycluster)
  mycluster <- gsub(".loom","",mycluster)
  outpath <- paste(WD,"/gene_passExpFilter_15cluster_",mycluster,".tsv", sep="")
  
  gene2keep %>%
    filter(cluster==mycluster) %>%
    fwrite(outpath, col.names=T, sep="\t")
}

#--- end ---#