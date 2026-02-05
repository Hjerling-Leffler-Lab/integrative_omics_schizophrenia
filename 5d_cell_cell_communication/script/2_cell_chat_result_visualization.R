#about: cell-cell communication network visualization with cellChat v2, Table T6: T6_sign_ligands_and_receptors.xlsx
#author: Lisa Bast
#date: 11.03.25
#version: 0.0.1

library(CellChat)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(future)
library(gridExtra)
library(openxlsx)

options(stringsAsFactors = FALSE)

opt_server <-  TRUE
opt_data_base = "original" #"manually curated"

#define pathes and load functions
code_path = getwd()
source("utils.R")

setwd("../")
main_folder = getwd()

results_path <- paste0(main_folder,"/output/")

## load objects
cellchat_SCZ = readRDS(paste0(results_path,"cellchat_SCZ_after_filtering_DB_",opt_data_base,".rds"))
cellchat_CTRL = readRDS(paste0(results_path,"cellchat_CTRL_after_filtering_DB_",opt_data_base,".rds"))

#compute network centrality
cellchat_SCZ <- netAnalysis_computeCentrality(cellchat_SCZ)
cellchat_CTRL <- netAnalysis_computeCentrality(cellchat_CTRL)
cellchat_SCZ <- computeCommunProbPathway(cellchat_SCZ)
cellchat_CTRL <- computeCommunProbPathway(cellchat_CTRL)

#merge objects:
object.list <- list(CTRL = cellchat_CTRL,SCZ = cellchat_SCZ)
cellchat_combined <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix=TRUE)

## plots
#ingoing vs outgoing interaction strength per cell type
visualize_ingoing_vs_outgoing_interaction_strength_per_celltype(cellchat_SCZ,cellchat_CTRL,results_path)
#to do overlay


#differential signaling: DEG analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "SCZ"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
#thresh.pc: 0.1 Threshold of the percent of cells expressed in one cluster
#thresh.fc: Threshold of Log Fold Change
#thresh.p: 0.05 Threshold of p-values
cellchat_combined <- identifyOverExpressedGenes(cellchat_combined, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.08, thresh.fc = 0.1, thresh.p = 0.05)

overexpressed_genes_sign_res <- cellchat_combined@var.features[[paste0(features.name, ".info")]][cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues.adj<0.05,]

#plot up-and down-reg signaling in schizophrenia
#differential interaction strength and number of interactions
visualize_differential_number_of_interactions_and_interaction_strength(cellchat_combined, results_path,"pathways")

#use the adjusted p-values for the selection:
cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues <- cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues.adj

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_combined, features.name = features.name)#, thresh=0.05)


# extract the ligand-receptor pairs with sign ligands in SCZ
net.ligand <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ",ligand.pvalues = 0.05)#,ligand.logFC = 0.15)#,thresh = 0.05)
# extract the ligand-receptor pairs with sign receptors in SCZ
net.receptor <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ",receptor.pvalues = 0.05)#,receptor.logFC = 0.15)#,thresh = 0.05)

browser()
save_sign_ligand_and_target_interactions_T6(net.ligand,net.receptor,results_path)


##apply different thresholds to visualize only sign communication with most extreme log2FCs
#sort data frames by logFC
#extract unique interaction_name values from top 20 rows
opt_filter_top5 <- TRUE
visualize_DEGS_chord_gene(cellchat_combined,object.list,results_path,net,net.ligand,"ligand",opt_filter_top5)
visualize_DEGS_chord_gene(cellchat_combined,object.list,results_path,net,net.receptor,"receptor",opt_filter_top5)


