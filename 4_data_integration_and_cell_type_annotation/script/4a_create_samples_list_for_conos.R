#Clustering snRNAseq - Human Cortex SCZ - Conos
#Create List with Donor expression
#Author: José A. Martínez-López
#Created on April 27, 2021
#Updated Dec 17, 2021

library(loomR)
library(Matrix)

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()

data_path <- paste0(main_project_path,'/3_quality_control/output/filtered/')
results_path <- paste0(main_path,'/output/Conos_clustering/metadata/')

#Input Loom file
file<-"Samples_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated.loom"
#number of samples:
N_SAMPLES<-85

#Connect Loom File
data_HC <- connect(paste0(data_path,file),mode = "r+", skip.validate = TRUE)

genes<-data_HC[['row_attrs/Gene']][]
Accession<-data_HC[['row_attrs/Accession']][]

genes_unique<-make.unique(genes)

CellID<-data_HC[['col_attrs/CellID']][]
Disease<-data_HC[['col_attrs/Disease']][]
Donor<-data_HC[['col_attrs/Donor']][]
Age<-data_HC[['col_attrs/Age']][]
scmap_ann_51CT<-data_HC[['col_attrs/CT_ann_ABM_MCA_scmap_cell2cluster_51_CTs']][]
scmap_ann_76CT<-data_HC[['col_attrs/CT_ann_ABM_MCA_scmap_cell2cluster_76_CTs']][]
scmap_score_51CT<-data_HC[['col_attrs/CT_ann_score_ABM_MCA_scmap_cell2cluster_51_CTs']][]
scmap_score_76CT<-data_HC[['col_attrs/CT_ann_score_ABM_MCA_scmap_cell2cluster_76_CTs']][]

names(Disease)<-CellID
names(Donor)<-CellID
names(Age)<-CellID
names(scmap_ann_51CT)<-CellID
names(scmap_ann_76CT)<-CellID
names(scmap_score_51CT)<-CellID
names(scmap_score_76CT)<-CellID


Disease<-factor(Disease)
Donor<-factor(Donor)
Age<-factor(Age)
scmap_ann_51CT<-factor(scmap_ann_51CT)
scmap_ann_76CT<-factor(scmap_ann_76CT)
scmap_score_51CT<-factor(scmap_score_51CT)
scmap_score_76CT<-factor(scmap_score_76CT)

samples<-list()

for (i in 1:N_SAMPLES) {
    
    s<-data_HC[['col_attrs/Donor']][]==paste0("S",i)
    
    if(sum(s)>0){    
        CellIDsample<-data_HC[['col_attrs/CellID']][s]
        S<-t(data_HC[['matrix']][s,])
        row.names(S)<-genes
        colnames(S)<-CellIDsample
        samples[[paste0("S",i)]]<-as(S, "dgCMatrix")
    }
}

saveRDS(genes,file=paste0(results_path,"Genes.rds"))
saveRDS(Accession,file=paste0(results_path,"Accession.rds"))
saveRDS(samples, file =paste0(results_path, "SamplesList_TH_and_D_adj_filtered_and_CT_annotated_20211217.rds"))
saveRDS(Disease,file = paste0(results_path,"Disease.rds"))
saveRDS(Donor,file = paste0(results_path,"Donor.rds"))
saveRDS(Age,file = paste0(results_path,"Age.rds"))
saveRDS(scmap_ann_51CT,file = "scmap_ann_51CT.rds")
saveRDS(scmap_ann_76CT,file = "scmap_ann_76CT.rds")
saveRDS(scmap_score_51CT,file = "scmap_score_51CT.rds")
saveRDS(scmap_score_76CT,file = "scmap_score_76CT.rds")


