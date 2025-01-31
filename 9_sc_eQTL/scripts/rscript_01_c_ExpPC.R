#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(tidyverse)
library(FactoMineR)
library(psych)

# read expression data
mylvl=args[1]
WD <- "./9_sc_eQTL/data/covar/" # specify
setwd(WD)
myfiles=list.files(WD, pattern="_agg_gene_donor_count.tsv")

# loop over CTs
for(i in 1:length(myfiles)){
  myCT=gsub("_agg_gene_donor_count.tsv","",myfiles[i])
  d <- fread(paste(WD,"/",myCT,"_agg_gene_donor_count.tsv",sep=""))
  e <- d %>%
    group_by(Donor) %>%
    mutate(donor_sum=sum(agg_count),
           tpm=agg_count*1000000/donor_sum) %>%
    ungroup() %>%
    select(-donor_sum)
  
  f1 <- e %>% pivot_wider(id_cols = Donor, names_from = ENSGID, values_from = tpm)
  #- PCA 
  res.pca1 <- PCA(t(f1[,-1]),  graph = F)
  eigenvalues <- res.pca1$eig %>% as.data.frame()
  eigenvalues$PC <- gsub("comp ","",rownames(eigenvalues)) %>% as.integer()
  
  if(i==1){df <- eigenvalues %>% mutate(CT=myCT)}else{
    df <- rbind(df,
                eigenvalues %>% mutate(CT=myCT))
  }
}

fwrite(df,paste("df.Exp.PC_",mylvl,".tsv",sep="")
       ,sep="\t", col.names = T)
