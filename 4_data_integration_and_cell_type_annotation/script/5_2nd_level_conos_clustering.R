#about: Clustering snRNAseq - Conos - 2nd level: excitatory , inhibitory and non-neuronal cells seperately, Fig. S3A
#Author: José A. Martínez-López, restructured and merged by Lisa Bast
#date: Created on Dec 6, 2020 updated on Dec 9,2021
#version: 0.0.2

library(pagoda2)
library(conos)
library(Matrix)

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()
metadata_path <- paste0(main_path,"/output/Conos_clustering/metadata/")
NCORES<-12

Disease<-readRDS(file=paste0(metadata_path,"Disease.rds"))
Donor<-readRDS(file=paste0(metadata_path,"Donor.rds"))
Age<-readRDS(file=paste0(metadata_path,"Age.rds"))
scmap_ann<-readRDS(file=paste0(metadata_path,"scmap_ann_76CT.rds"))
scmap_score<-readRDS(file=paste0(metadata_path,"scmap_score_76CT.rds"))

for (ct_class in c("excitatory","inhibitory","nonneuronal")){
  if (ct_class=="excitatory"){
    subfolder="Excitatory"
  } else if (ct_class=="inhibitory"){
    subfolder="Inhibitory"
  } else if (ct_class=="nonneuronal"){
    subfolder="Nonneuronal"
  }
  data_path <- paste0(main_path,"/output/Conos_clustering/1st_level/")
  results_path <- paste0(main_path,'/output/Conos_clustering/2nd_level/',subfolder,"/")

  samples<-readRDS(file=paste0(data_path,"samples_",ct_class,"_cells_20211218.rds"))
  
  str(samples, 1)
  
  #Preprocessing of samples
  panel.preprocessed <- lapply(samples, basicP2proc, n.cores=NCORES, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
  
  #Create Conos object
  con <- Conos$new(panel.preprocessed, n.cores=NCORES)
  
  #Remove original data and Pagoda object
  rm(samples,panel.preprocessed)
  
  str(con$samples,1)
  
  #Clusters of each sample
  plot01_samples<-con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)
  
  #Build Graph
  con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=FALSE, verbose=TRUE)
  
  #Global clusters
  con$findCommunities(method=leiden.community, resolution=2)
  
  plotA<-con$plotGraph(alpha=0.1)
  plotB<-con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
  plotC<-con$plotGraph(mark.groups=TRUE, alpha=0.1, show.legend=TRUE)
  plotD<-con$plotGraph(gene='SLC17A7')
  plotE<-con$plotGraph(gene='GAD1')
  
  saveRDS(plotA,file=paste0(results_path,"plotA.rds"))
  saveRDS(plotB,file=paste0(results_path,"plotB.rds"))
  saveRDS(plotC,file=paste0(results_path,"plotC.rds"))
  saveRDS(plotD,file=paste0(results_path,"plotD.rds"))
  saveRDS(plotE,file=paste0(results_path,"plotE.rds"))
  
  #Clusters of each sample (same colors)
  plot02_samplesSameClusters<-con$plotPanel(font.size=4)
  
  #Composition of clusters
  plot03_composition<-plotClusterBarplots(con, legend.height = 0.1)
  
  saveRDS(plot01_samples,file=paste0(results_path,"plot01.rds"))
  saveRDS(plot02_samplesSameClusters,file=paste0(results_path,"plot02.rds"))
  saveRDS(plot03_composition,file=paste0(results_path,"plot03.rds"))
  #saveConosForScanPy(con, output.path=outputpath, hdf5_filename="CONOS_analysis01.h5", verbose=TRUE)
  
  #Parameters Graph
  con$embedGraph(method="UMAP",n.cores=NCORES)
  
  #Global Graph
  plot04_global<-con$plotGraph(size=0.1)
  saveRDS(plot04_global,file=paste0(results_path,"plot04.rds"))
  pdf(paste0(results_path,"plot04_global",".pdf"),width=20,height=20)
  plot(plot04_global)
  dev.off()
  
  plotSamples<-con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
  saveRDS(plotSamples,file=paste0(results_path,"plotSamples.rds"))
  pdf(paste0(results_path,"plotSamples",".pdf"),width=20,height=20)
  plot(plotSamples)
  dev.off()
  
  plot_disease<-con$plotGraph(groups=Disease,mark.groups=FALSE,alpha=0.2,plot.na=FALSE,title='Disease',show.legend=TRUE)
  saveRDS(plot_disease,file=paste0(results_path,"plot_disease.rds"))
  pdf(paste0(results_path,"plot_disease",".pdf"),width=20,height=20)
  plot(plot_disease)
  dev.off()
  
  plot_donor<-con$plotGraph(groups=Donor,mark.groups=FALSE,alpha=0.2,plot.na=FALSE,title='Donor',show.legend=TRUE)
  saveRDS(plot_donor,file=paste0(results_path,"plot_donor.rds"))
  pdf(paste0(results_path,"plot_donor",".pdf"),width=20,height=20)
  plot(plot_donor)
  dev.off()
  
  
  plot_age<-con$plotGraph(groups=Age,mark.groups=FALSE,alpha=0.2,plot.na=FALSE,title='Age',show.legend=TRUE)
  saveRDS(plot_age,file=paste0(results_path,"plot_age.rds"))
  pdf(paste0(results_path,"plot_age",".pdf"),width=20,height=20)
  plot(plot_age)
  dev.off()
  
  plot_scmap_ann<-con$plotGraph(groups=scmap_ann,mark.groups=TRUE,alpha=0.2,plot.na=FALSE,title='scmap_ann',show.legend=FALSE)
  saveRDS(plot_scmap_ann,file=paste0(results_path,"plot_scmap_ann.rds"))
  pdf(paste0(results_path,"plot_scmap_ann",".pdf"),width=20,height=20)
  plot(plot_scmap_ann)
  dev.off()
  
  plot_scmap_score<-con$plotGraph(groups=scmap_score,mark.groups=FALSE,alpha=0.2,plot.na=TRUE,title='scmap_score',show.legend=TRUE)
  saveRDS(plot_scmap_score,file=paste0(results_path,"plot_scmap_score.rds"))
  #pdf(paste0("plot_scmap_score",".pdf"),width=20,height=20)
  #plot(plot_scmap_score)
  #dev.off()
  
  #saveConosForScanPy(con, output.path=outputpath, hdf5_filename="CONOS_analysis_leiden_resolution5.h5", verbose=TRUE)
  
  metadata <- data.frame(Cluster=con$clusters$leiden$groups)
  saveRDS(metadata,file=paste0(results_path,"clustersAnnotation.rds"))
  
  #markers from literature
  plot05<-con$plotGraph(gene='SLC17A7')
  saveRDS(plot05,file=paste0(results_path,"plot05.rds"))
  
  plot06<-con$plotGraph(gene='GAD1')
  saveRDS(plot06,file=paste0(results_path,"plot06.rds"))
  
  plot07<-con$plotGraph(gene='LAMP5')
  saveRDS(plot07,file=paste0(results_path,"plot07.rds"))
  
  plot08<-con$plotGraph(gene='PVALB')
  saveRDS(plot08,file=paste0(results_path,"plot08.rds"))
  
  plot09<-con$plotGraph(gene='VIP')
  saveRDS(plot09,file=paste0(results_path,"plot09.rds"))
  
  plot10<-con$plotGraph(gene='SST')
  saveRDS(plot10,file=paste0(results_path,"plot10.rds"))
  
  plot11<-con$plotGraph(gene='SV2C')
  saveRDS(plot11,file=paste0(results_path,"plot11.rds"))
  
  plot12<-con$plotGraph(gene='CUX2')
  saveRDS(plot12,file=paste0(results_path,"plot12.rds"))
  
  plot13<-con$plotGraph(gene='RORB')
  saveRDS(plot13,file=paste0(results_path,"plot13.rds"))
  
  plot14<-con$plotGraph(gene='AQP4')
  saveRDS(plot14,file=paste0(results_path,"plot14.rds"))
  
  plot15<-con$plotGraph(gene='MOG')
  saveRDS(plot15,file=paste0(results_path,"plot15.rds"))
  
  plot16<-con$plotGraph(gene='CSF1R')
  saveRDS(plot16,file=paste0(results_path,"plot16.rds"))
  
  plot17<-con$plotGraph(gene='PDGFRA')
  saveRDS(plot17,file=paste0(results_path,"plot17.rds"))
  
  #markers
  #de.info <- con$getDifferentialGenes(groups=con$clusters$leiden$groups,append.auc=TRUE)
  #saveRDS(de.info,file="de.rds")
  
  #heatmapplot<-plotDEheatmap(con,con$clusters$leiden$groups,de.info,n.genes.per.cluster=5,row.label.font.size=7)
  #saveRDS(heatmapplot,file="heatm#ap_plot.rds")
  
}

