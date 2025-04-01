#date: 24.03.2021
#author: Lisa Bast
#version: 1.0.1
#plot heatmap with seurat
#R version: [Default] [64-bit] C:\Program Files\R\R-4.0.2

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()
results_folder_str <- paste0(main_path,"/output/CT_anno_scmap/ABM_MCA/")
add_str_pagoda = '_pagoda' #''

library(loomR)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(patchwork)
#library(Matrix)
library(iotools) # required for read.csv.raw
library(ggplot2) #required for figure cosmetics
library(stringr)

#n_CT=76 # 51

#Input Loom file
dir<-paste0(main_path,'/output/CT_specific_files/3_classes/')
files<-c(paste0("Samples_NonNeuronal_conos",add_str_pagoda,"_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"),
         paste0("Samples_Excitatory_conos",add_str_pagoda,"_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"),
         paste0("Samples_Inhibitory_conos",add_str_pagoda,"_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"),
         paste0("Samples",add_str_pagoda,"_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"))

for (file in files){
	filename_pieces <- str_split(file[[1]], '_')
	#Connect Loom File
	data_HC <- connect(paste0(dir,file),mode = "r+", skip.validate = TRUE)
	#S_object <- as.Seurat(data_HC,cells="CellID",features="Gene")

	genes<-data_HC[['row_attrs/Gene']][]
	genes<-make.unique(genes)

	CellID<-data_HC[['col_attrs/CellID']][]
	#
	#Donor<-data_HC[['col_attrs/Donor']][]
	#Library<-data_HC[['col_attrs/Library']][]
	#PMI_h<-data_HC[['col_attrs/PMI_h']][]
	#Sex<-data_HC[['col_attrs/Sex']][]
	#mean_counts<-data_HC[['col_attrs/mean_counts_per_barcode']][]
	#Age<-data_HC[['col_attrs/Age']][]
	#scmap_ann<-data_HC[[paste0(paste0('col_attrs/CT_ann_ABM_MCA_scmap_cluster_',n_CT),'CTs')]][]
	#scmap_score<-data_HC[[paste0(paste0('col_attrs/CT_ann_score_ABM_MCA_scmap_cluster_',n_CT),'CTs')]][]
	if ((filename_pieces[[1]][2]=='Inhibitory') || (filename_pieces[[1]][2]=='Excitatory') || (filename_pieces[[1]][2]=='NonNeuronal') || (filename_pieces[[1]][1]=="S1")){
    conos_cluster <- as.numeric(data_HC[['col_attrs/Conos_2nd_level']][])
	} else{
	  conos_cluster <- as.numeric(data_HC[['col_attrs/Conos_2nd_level']][])
	}

	full.matrix <- t(data_HC[["matrix"]][, ])
	colnames(full.matrix) = CellID
	rownames(full.matrix) = genes

	#genes used for CT annotation with scmap:
	#gene_list_scmap_path <- paste0(main_path,"data/reference_data_sets/ABM_MCA/red_to_2000_HVG_with_cell_ranger/")
	#setwd(gene_list_scmap_path)
	#matrix_scmap <- read.csv.raw(file="9606_symbol.csv")
	#gene_list_scmap <- matrix_scmap$symbol

	#gene_list_intersect <- intersect(rownames(full.matrix),gene_list_scmap)
	#reduced.matrix <- full.matrix[gene_list_intersect,]

	S_object <- CreateSeuratObject(counts = full.matrix, project = "SCZ_human_PFC", assay="RNA", min.cells = 0, min.features = 0, names.field = 1)
	S_object <- AddMetaData(S_object,conos_cluster,col.name ='conos_cluster')
	#S_object <- AddMetaData(S_object,Disease,col.name = 'Disease')
	#S_object <- AddMetaData(S_object,Donor,col.name = 'Donor')
	#S_object <- AddMetaData(S_object,Library,col.name = 'Library')
	#S_object <- AddMetaData(S_object,PMI_h,col.name = 'PMI_h')
	#S_object <- AddMetaData(S_object,Sex,col.name = 'Sex')
	#S_object <- AddMetaData(S_object,Age,col.name = 'Age')
	#S_object <- AddMetaData(S_object,scmap_ann,col.name = 'scmap_ann')
	#S_object <- AddMetaData(S_object,scmap_score,col.name = 'scmap_score')

	S_object <- NormalizeData(S_object, normalization.method = "LogNormalize", scale.factor = 10000)

	#determine highly variable genes:
	S_object <- FindVariableFeatures(S_object, selection.method = "vst", nfeatures = 3000)
	#VariableFeatures(S_object)
	all.genes <- rownames(S_object)

	#scale data taking into account all genes:
	S_object <- ScaleData(S_object, features = all.genes)
	dim(S_object@assays$RNA@scale.data)

	#perform PCA on highly variable genes
	#S_object  <- RunPCA(S_object , features = VariableFeatures(object = S_object))
	#S_object <- FindNeighbors(S_object, dims = 1:10)
	#S_object <- FindClusters(S_object, resolution = 0.5)
	#S_object.markers <- FindAllMarkers(S_object, test.use='wilcox')

	#top_20 <- S_object.markers %>% group_by(conos_cluster) %>% top_n(n = 20, wt = avg_log2FC)

	#determine for each conos cluser the most differentially expressed genes (compared to the other clusters) and store all significant genes with their log2 FC

	#for each conos cluster specific gene list:
	#identify number of genes proportional to cells in that cluster 
	#take genes with the highest absolute logFC values 
	#order genes within a cluster by logFC values (neg to pos)
	n_target_genes <-3000
	n_cells <- length(conos_cluster)
	for (i in sort(unique(conos_cluster))){
	  print(i)
	  n_cells_per_cluster <-sum(conos_cluster==i)
	  n_genes_per_cluster <- round(n_cells_per_cluster/n_cells*n_target_genes)
	  marker_genes_i <- FindMarkers(S_object,CellID[conos_cluster==i],test.use='wilcox')
	  marker_genes_i$conos_cluster <- c(rep(i,length(rownames(marker_genes_i))))
	  #add column with absolute logFC values
	  marker_genes_i$abs_avg_log2FC <- abs(marker_genes_i$avg_log2FC)
	  #sort dataframe according to absolute logFC values and extract the best (TH_abs_avg_log2FC_i)th genes
	  marker_genes_i_sorted_abs <- marker_genes_i[order(marker_genes_i$abs_avg_log2FC,decreasing = TRUE),]
	  gene_list_per_cluster <- rownames(marker_genes_i_sorted_abs[1:n_genes_per_cluster,])
	  #add new column to dataframe stating if gene is selected or not
	  marker_genes_i$selected <- rep(FALSE,dim(marker_genes_i)[1])
	  marker_genes_i[gene_list_per_cluster,]$selected<-TRUE
	  #order dataframe according to log FC
	  marker_genes_i_sorted <- marker_genes_i[order(marker_genes_i$avg_log2FC,decreasing = FALSE),]
	  sorted_gene_list_per_cluster <- rownames(marker_genes_i_sorted[marker_genes_i_sorted$selected==TRUE,])
	  n_genes_per_cluster_downregulated <- length(rownames(marker_genes_i_sorted[marker_genes_i_sorted$selected==TRUE & marker_genes_i_sorted$avg_log2FC<0,]))
	  n_genes_per_cluster_upregulated <- length(rownames(marker_genes_i_sorted[marker_genes_i_sorted$selected==TRUE & marker_genes_i_sorted$avg_log2FC>0,]))
	  if (i==1){
		  marker_genes_per_conos_cluster <- marker_genes_i_sorted
		  marker_genes_list <- sorted_gene_list_per_cluster
		  marker_genes_list_info <- c(rep(paste0(paste0("Conos_cluster_",as.character(i)),"_downregulated"),n_genes_per_cluster_downregulated),rep(paste0(paste0("Conos_cluster_",as.character(i)),"_upregulated"),n_genes_per_cluster_upregulated))
	  } else{
		  marker_genes_per_conos_cluster <- rbind(marker_genes_per_conos_cluster,marker_genes_i_sorted,sort=FALSE)
		  marker_genes_list <- c(marker_genes_list,sorted_gene_list_per_cluster)
		  marker_genes_list_info <- c(marker_genes_list_info,rep(paste0(paste0("Conos_cluster_",as.character(i)),"_downregulated"),n_genes_per_cluster_downregulated),rep(paste0(paste0("Conos_cluster_",as.character(i)),"_upregulated"),n_genes_per_cluster_upregulated))
	  }
	}
	gene_list.df <- data.frame(marker_genes_list,marker_genes_list_info)
	
	
	#save list
	write.csv(gene_list.df, paste0(results_folder_str,"/gene_list_",filename_pieces[[1]][2],".csv"), row.names = TRUE)
	
	#TH_abs_avg_log2FC_i <- sort(abs(marker_genes_per_conos_cluster$avg_log2FC[marker_genes_per_conos_cluster$conos_cluster==1]),decreasing = TRUE)[n_genes_per_cluster]
	graphic_filename_1 <- paste0(results_folder_str,"/heatmap_unscaled_",filename_pieces[[1]][1],'_',filename_pieces[[1]][2],".pdf")
	pdf(file = graphic_filename_1,   # The directory you want to save the file in
		width = 6, # The width of the plot in inches
		height = 6) # The height of the plot in inches
	print(DoHeatmap(object=S_object, features=marker_genes_list, group.by="conos_cluster", slot='data', label=TRUE, size=1, angle=90)+
	  theme(axis.text.x = element_text(size = 3)))
	dev.off()
	
	graphic_filename_2 <- paste0(results_folder_str,"/heatmap_scaled_",filename_pieces[[1]][1],'_',filename_pieces[[1]][2],".pdf")
	pdf(file = graphic_filename_2,   # The directory you want to save the file in
		width = 6, # The width of the plot in inches
		height = 6) # The height of the plot in inches
	print(DoHeatmap(object=S_object, features=marker_genes_list, group.by="conos_cluster", slot='scale.data', label=TRUE, size=1, angle=90)+
	  theme(axis.text.x = element_text(size = 3)))
	dev.off()
}

