#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# File name: rscript_01_b2_Get_Exp_Input_3CT.R
# Aim: generate expression input for eQTL analysis.

#--- script starts ---#

# Libraries
library(tidyverse)
library(data.table)
library(rhdf5)
library(snow)

#- specify paths
filepath <- "" # specify loom files (15CT)
filelist0 <- grep(list.files(path=filepath, pattern='\\.loom$'), pattern='downsampled', invert=TRUE, value=TRUE)
WD <- "./9_sc_eQTL/data/eQTL/pheno/3CT" # specify working directory
setwd(WD)

#- specify cell type
myCT=args[1]

#- read in geneMatrix 
# (with hg19 gene coordinates; filter for protein-coding genes on autosomes)
geneMx <- fread("geneMatrix.tsv.gz",stringsAsFactors=F,data.table=F) %>%
  filter(gene_type=="protein_coding",
         is.na(hg19g0)==FALSE,
         is.na(g1)==FALSE,
         is.na(g2)==FALSE,
         PAR=="FALSE",
         !hg19g0%in%c("chrM","chrX","chrY")) %>%
  mutate(tss=as.integer(ifelse(gstr=="+", g1, g2)),
         tss0=as.integer(tss-1)) %>%
  select(ensgid, hg19g0, tss0, tss, gstr)

#- read sample information (ID and batch)
sample.df <- fread("./workflow/brn_GRS_EUR_donorIDinternal8.tsv", # get through .R script ./9_sc_eQTL/scripts/99_brn_GRS_EUR_donorIDinternal8.R
                   stringsAsFactors = F, data.table=F) %>%
  select(donor_ID_internal_8, IID)

#- 1. Get coordinates for QC:ed genes
genelist <- fread(paste(WD,"/gene_passExpFilter_3cluster_",myCT,".tsv", sep=""),
                  stringsAsFactors=F, data.table=F)

geneMx <- geneMx %>%
  filter(ensgid %in% genelist$ENSGID)

#- 2. Get expression for QC:ed genes
mypattern=paste("_",myCT,"_",sep="")
filelist <- grep(list.files(filepath, pattern = mypattern), pattern='downsampled', invert=TRUE, value=TRUE)

for (i in 1:length(filelist)){
  file=paste(filepath, filelist[i], sep="")
  mycluster=basename(file)
  mycluster <- gsub("Samples_cellranger_pagoda_TH_and_D_adj_filtered_conos_cluster_based_","",mycluster)
  mycluster <- gsub(".loom","",mycluster)
  
  h5f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
  cell_donor <- cbind(h5f$col_attrs$CellID,
                      h5f$col_attrs$Donor) %>%
    as.data.frame()
  colnames(cell_donor) <- c("CellID", "Donor")
  
  exp <- as.data.frame(t(h5f$matrix))
  colnames(exp) <- h5f$col_attrs$CellID
  exp$Accession <- h5f$row_attrs$Accession
  exp <- exp %>% 
    filter(Accession %in% geneMx$ensgid) %>%
    select(Accession, everything())
  exp.long0 <- exp %>% 
    pivot_longer(col=-Accession, 
                 names_to = "CellID",
                 values_to = "count")
  rm(exp)
  exp.long <- exp.long0 %>% 
    left_join(cell_donor, by="CellID") %>% 
    group_by(Accession, Donor) %>%
    summarise(agg_count=sum(count)) %>%
    ungroup()
  rm(exp.long0)
  if(i==1){
    my.exp.long <- exp.long}else{
      my.exp.long <- rbind(my.exp.long, exp.long) %>%
        group_by(Accession, Donor) %>%
        summarise(agg_count=sum(agg_count)) %>%
        ungroup()
    }
}

colnames(my.exp.long)[1] <- "ENSGID"

# output for calculating expression PCs (covariates of QTL analysis, see 01_c_prepare_covariates.Rmd)
# system(paste("mkdir ",WD,"/3CT/", sep=""))
#outfile <- paste(WD,"/3CT/", myCT, "_agg_gene_donor_count.tsv", sep="")
#fwrite(my.exp.long, outfile, sep="\t", col.name=T)

# convert to tmp
my.exp.long <- my.exp.long %>%
  group_by(Donor) %>%
  mutate(total_count=sum(agg_count))%>%
  mutate(exp_tpm=agg_count*1000000/total_count) %>%
  ungroup()

#- 3. match to the correct ID and reorder
my.exp.long <- my.exp.long %>%
  left_join(sample.df, by=c("Donor"="donor_ID_internal_8")) %>%
  filter(is.na(IID)==F, IID!="") %>% # keep only samples of EUR ancestry
  mutate(batch=ifelse(substr(IID,4,4)=="-","batch1", "batch2")) %>%
  arrange(batch, IID, ENSGID)

df <- my.exp.long %>%
  select(ENSGID, IID, exp_tpm) %>% # filter for useful columns
  pivot_wider(names_from = IID, 
              values_from = exp_tpm) %>% # to wide format
  inner_join(geneMx, by=c("ENSGID"="ensgid")) %>% # map coordinates
  mutate(ENSGID1=ENSGID) %>%
  select(hg19g0, tss0, tss, ENSGID1, ENSGID, gstr, everything()) %>%
  arrange(hg19g0, tss0, tss)

colnames(df)[1:6] <- c("#Chr", "start", "end", "pid", "gid", "strand")

#- output per chr
for (i in 1:22){
  dfi <- df %>%
    filter(`#Chr`==paste0("chr",i))
  fwrite(dfi,
         paste0(myCT,"_exp_tpm.chr",i,"_TSS.bed"),
         sep="\t",col.names = T)
}

#--- end ---#