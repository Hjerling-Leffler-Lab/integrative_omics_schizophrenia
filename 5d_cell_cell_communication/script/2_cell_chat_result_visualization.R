#about: cell-cell communication network visualization with cellChat v2
#author: Lisa Bast
#date: 11.03.25
#version: 0.0.1

library(CellChat)
library(patchwork)
library(stringr)
library(RColorBrewer)
library(future)
library(gridExtra)

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


## visualize results
#aggregated_network
visualize_aggr_network(cellchat_combined,cellchat_SCZ,cellchat_CTRL,results_path)

#compare total number of interactions:
visualize_number_of_interactions_and_interaction_strength(cellchat_combined, results_path)


#differential number of interactions/ interaction strength
visualize_differential_number_of_interactions_and_interaction_strength(cellchat_combined, results_path)

#ingoing vs outgoing interaction strength per cell type
visualize_ingoing_vs_outgoing_interaction_strength_per_celltype(cellchat_SCZ,cellchat_CTRL,results_path)

#Identify signaling changes associated with one cell group
visualize_ingoing_vs_outgoing_signaling_changes_per_celltype(cellchat_combined,results_path)

#for interesting pairs of cell types:
#CTs <- levels(cellchat_combined@idents$joint)
# [1] "Astrocytes"         "Endothelial/mural"  "Exc L2-3 IT I"      "Exc L2-3 IT II"     "Exc L3-4 IT"        "Exc L3-6 IT"        "Exc L5-6 CT and NP" "Exc L5-6 IT I"     
#[9] "Exc L5-6 IT II"     "Inh LAMP5"          "Inh PVALB"          "Inh SST"            "Inh VIP"            "Microglia"          "Oligodendrocytes"   "OPCs"    
source_celltypes <- c("Exc L5-6 IT II","Inh LAMP5",#weaker interaction
                      "Exc L2-3 IT I","Exc L3-6 IT","Exc L5-6 IT II" ,"Exc L5-6 IT II",#stronger interaction
                      "Endothelial/mural","Exc L3-4 IT",#more interactions
                      "OPCs","OPCs","OPCs","Exc L2-3 IT I" ,"Exc L2-3 IT I" ,"Exc L2-3 IT I" ,"Exc L2-3 IT I" ,"Exc L2-3 IT I")#fewer interactions
target_celltypes <- c("Oligodendrocytes","Oligodendrocytes",#weaker interaction
                      "Inh LAMP5","Inh LAMP5","Exc L2-3 IT I","Exc L3-6 IT",#stronger interaction
                      "Astrocytes","Astrocytes",#more interactions
                      "Exc L2-3 IT I","Exc L3-4 IT","Inh LAMP5","Endothelial/mural","Exc L2-3 IT I","Exc L5-6 IT II","Inh LAMP5","OPCs")#fewer interactions
#ranked signaling pathways
visualize_ranked_signaling_pathways_between_CT_pairs(cellchat_combined,results_path,source_celltypes,target_celltypes)
  
#dysfunctional signaling: comparison of communication probabilities
for (i in seq(1,length(source_celltypes))){
  visualize_difference_in_communication_probabilities(cellchat_combined,source_celltypes[i],target_celltypes[i],results_path)
}

#dysfunctional signaling: DEG analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "SCZ"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
#thresh.pc: 0.1 Threshold of the percent of cells expressed in one cluster
#thresh.fc: Threshold of Log Fold Change
#thresh.p: 0.05 Threshold of p-values
cellchat_combined <- identifyOverExpressedGenes(cellchat_combined, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.2, thresh.p = 0.05)

overexpressed_genes_sign_res <- cellchat_combined@var.features[[paste0(features.name, ".info")]][cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues.adj<0.05,]
write.csv(overexpressed_genes_sign_res,paste0(results_path,"overexpressed_genes.csv"), row.names = FALSE)
  
#use the adjusted p-values for the selection:
cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues <- cellchat_combined@var.features[[paste0(features.name, ".info")]]$pvalues.adj

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_combined, features.name = features.name, thresh=0.05)

# extract the ligand-receptor pairs with upregulated ligands in SCZ
net.up <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ",ligand.pvalues = 0.05,receptor.pvalues = 0.05,thresh = 0.05)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in CTRL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_combined, net = net, datasets = "CTRL",ligand.pvalues = 0.05, receptor.pvalues = 0.05,thresh = 0.05)

write.csv(net.up, paste0(results_path,"down_upgulated_interactions.csv"), row.names = FALSE)
write.csv(net.down, paste0(results_path,"down_regulated_interactions.csv"), row.names = FALSE)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_combined)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_combined)

pathways_all <- get_pathways(cellchat_SCZ,cellchat_CTRL)
for (i in seq(1,length(source_celltypes))){
  visualize_DEGS_bubble(cellchat_combined,object.list,source_celltypes[i],target_celltypes[i],results_path)
}

visualize_DEGS_chord_gene(cellchat_combined,object.list,results_path,net.up,net.down)

pathways.up <- unique(net.up$pathway_name)
pathways.down <- unique(net.down$pathway_name)
for (pw in pathways.up){
  visualize_DEGS_chord_cell(cellchat_combined,object.list,results_path,net.up,net.down,pw,"up")
}
for (pw in pathways.down){
  visualize_DEGS_chord_cell(cellchat_combined,object.list,results_path,net.up,net.down,pw,"down")
}

#visualize_DEGS(cellchat_combined,object.list,source_celltypes[i],target_celltypes[i],results_path,"bubble",net.up,net.down,pathways_all)
#visualize_DEGS(cellchat_combined,object.list,source_celltypes[i],target_celltypes[i],results_path,"chord",net.up,net.down,pathways_all)


#network for specific pathways
pathways.show <- c("NRXN") #"SOMATOSTATIN", "VIP" 
visualize_cell_cell_comm_network(cellchat_SCZ,results_path,gene.up,"SCZ")
visualize_cell_cell_comm_network(cellchat_CTRL,results_path,pathways.show,"CTRL")

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human')

# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human')

#Differential number of interactions or interaction strength among different cell types
#"Inh LAMP5" "Oligodendrocytes"  "OPCs" "Astrocytes" "Endothelial/mural"
#group.cellType <- c(rep("Exc L2-3 IT I", 4), rep("Exc L3-6 IT", 4), rep("Exc L5-6 IT II", 4))
#group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
#object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
#par(mfrow = c(1,2), xpd=TRUE)
#for (i in 1:length(object.list)) {
#  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
#}

#Identify signaling groups based on their functional similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "functional")
cellchat_combined <- netEmbedding(cellchat_combined, type = "functional")
#Error in runUMAP(Similarity, min_dist = min_dist, n_neighbors = n_neighbors,  : 
#Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).

cellchat_combined <- netClustering_multisession(cellchat_combined, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_combined, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat_combined, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat_combined <- computeNetSimilarityPairwise(cellchat_combined, type = "structural")
cellchat_combined <- netEmbedding(cellchat_combined, type = "structural")
cellchat_combined <- netClustering(cellchat_combined, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_combined, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_combined, type = "structural", nCol = 2)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat_combined, type = "functional")