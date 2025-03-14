#about:  Clustering snRNAseq - Human Cortex SCZ - Conos
#        for 2nd level clustering
#Author: José A. Martínez-López, restructured and merged by Lisa Bast
#date:   Created on Dec 9, 2021
#version:0.0.1

library(loomR)
library(Matrix)

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()


#Input Loom file
file<-"Samples_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated.loom"
N_SAMPLES<-85

#Connect Loom File
data_HC <- connect(paste0(data_path,file),mode = "r+", skip.validate = TRUE)

genes<-data_HC[['row_attrs/Gene']][]
#genes_unique<-make.unique(genes)

CellID<-data_HC[['col_attrs/CellID']][]

conos_clusters<-readRDS(file=paste0(first_level_results_path,"clustersAnnotation.rds"))

#Ordering clusters
conos_clusters_ordered<-data.matrix(conos_clusters[order(factor(row.names(conos_clusters),levels=CellID)),])

for (ct_class in c("excitatory","inhibitory","nonneuronal")){
  if (ct_class=="excitatory"){
    subfolder="Excitatory"
  } else if (ct_class=="inhibitory"){
    subfolder="Inhibitory"
  } else if (ct_class=="nonneuronal"){
    subfolder="Nonneuronal"
  }
  data_path <- paste0(main_project_path,'/3_quality_control/output/filtered/')
  results_path <- paste0(main_path,'/output/Conos_clustering/2nd_level/',subfolder,"/")
  first_level_results_path <- results_path
  NCORES<-12
  
  
  samples<-list()
  
  if (ct_class=="excitatory"){
    for (i in 1:N_SAMPLES) {
      s<-data_HC[['col_attrs/Donor']][]==paste0("S",i)
      c<-conos_clusters_ordered %in% c("2","8","7","17","43","13","6","20","18","24","39","4","27","9","3","14","37","21","32","26","34","48","35")
      
      selection<-s&c
      
      if (sum(s)>0){
        CellIDsample<-data_HC[['col_attrs/CellID']][selection]
        S<-t(data_HC[['matrix']][selection,])
        row.names(S)<-genes
        colnames(S)<-CellIDsample
        samples[[paste0("S",i)]]<-as(S, "dgCMatrix")
      }
    } 
  } else if (ct_class=="inhibitory"){
    for (i in 1:N_SAMPLES) {
      s<-data_HC[['col_attrs/Donor']][]==paste0("S",i)
      c<-conos_clusters_ordered %in% c("15","30","31","40","10","16","41","28","19","33","11","38","12","44","22","29","23","49")
      
      selection<-s&c
      
      if (sum(s)>0){
        CellIDsample<-data_HC[['col_attrs/CellID']][selection]
        S<-t(data_HC[['matrix']][selection,])
        row.names(S)<-genes
        colnames(S)<-CellIDsample
        samples[[paste0("S",i)]]<-as(S, "dgCMatrix")
      }
    }
    
  } else if (ct_class=="nonneuronal"){
    for (i in 1:N_SAMPLES) {
      s<-data_HC[['col_attrs/Donor']][]==paste0("S",i)
      c<-conos_clusters_ordered %in% c("1","25","47","45","5","42","36","46")
      
      selection<-s&c
      
      if (sum(s)>0){
        CellIDsample<-data_HC[['col_attrs/CellID']][selection]
        S<-t(data_HC[['matrix']][selection,])
        row.names(S)<-genes
        colnames(S)<-CellIDsample
        samples[[paste0("S",i)]]<-as(S, "dgCMatrix")
      }
    }
  }
  
  saveRDS(samples, file = paste0(results_path,"samples_",ct_class,"_cells_20211218.rds"))
  
  print("Done")
  
}
