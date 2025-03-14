# File name: 00_01_b1_Filter_Genes_3CTs.R
# Aim: Remove genes for eQTL analysis; keep genes with ≥1 count in ≥10% cells in each cell type

#--- script starts ---#

# libraries
library(tidyverse)
library(loomR)
library(data.table)

# specify paths
filepath <- "" # folder of loom files
filelist <- grep(list.files(path=filepath, pattern='\\.loom$'), pattern='downsampled', invert=TRUE, value=TRUE)
WD <- "./9_sc_eQTL/data/eQTL/pheno/3CT" # working directory
setwd(WD)

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

# aggregate per class
mydf$category <- gsub("\\_.*","",mydf$cluster)
unique(mydf$category)
mydf <- mydf %>%
  mutate(category=ifelse(category%in%c("Excitatory","Inhibitory"), category, "NonNeuronal"))
mydf2 <- mydf %>% 
  group_by(ENSGID, category) %>%
  summarise(count_cell_0exp=sum(count_cell_0exp),
            N_cell_total=sum(N_cell_total)) %>% 
  ungroup()

# genes pass expression filter
gene2keep <- mydf2 %>% filter(count_cell_0exp <= 0.9*N_cell_total)

# output
fwrite(gene2keep %>% filter(category=="Excitatory"), 
       paste(WD,"/gene_passExpFilter_3cluster_Excitatory.tsv", sep=""), 
       col.names = T, sep="\t")
fwrite(gene2keep %>% filter(category=="Inhibitory"), 
       paste(WD,"/gene_passExpFilter_3cluster_Inhibitory.tsv", sep=""), 
       col.names = T, sep="\t")
fwrite(gene2keep %>% filter(category=="NonNeuronal"), 
       paste(WD,"/gene_passExpFilter_3cluster_NonNeuronal.tsv", sep=""), 
       col.names = T, sep="\t")

#--- end ---#