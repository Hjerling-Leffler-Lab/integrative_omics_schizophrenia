#about: cell-cell communication inference with cellChat v2
#author: Lisa Bast
#date: 11.03.25
#version: 0.0.1


# to do_ make cell type names shorter: replace Inhibitory with Inh, Excitatory with Exc, Layer with L, neurons with " "
library(CellChat)
library(patchwork)
library(loomR)
library(future)
library(stringr)

options(stringsAsFactors = FALSE)

#better: paralellize for server
subset_conditions <- c("CTRL","SCZ")#,"both")
opt_data_bases = c("original", "manually curated")

opt_server <-  TRUE

if (opt_server==TRUE){
  args <- commandArgs(trailingOnly=TRUE)
}else{
  args <- "1"
}
n <- as.integer(args)
if (n%%2==0){
  i<-1
} else{
  i<-2
}
if (n<=2){
  j<-1
} else{
  j<-2
}
opt_data_base <- opt_data_bases[i]
subset_condition <- subset_conditions[j]

#define pathes and load functions
code_path = getwd()
source("utils.R")

setwd("../")
setwd("../")
main_path = getwd()
main_subfolder = "/5d_cell_cell_communication/"

#data_path<-paste0(main_path,'/4_data_integration_and_cell_type_annotation/output/')
if (opt_server==TRUE){
  data_path <- "/nas/depts/007/sullilab/projects/snrnaseq/data/filtered_loom_formatted_data_cellranger/"
  data_file <- "Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom" 
}else{
  data_path <- "D:/Documents/projects/SCZ_human_cortex/data/filtered_loom_formatted_data_cellranger/"
  data_file <- "Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered_subsampled_to_10_percent_cells.loom" 
}
results_path <- paste0(main_path,main_subfolder,"output/")


##1) load data and create cell_chat object:
cellchat <- load_data_and_get_cell_chat_object(data_path,data_file,subset_condition)
print("----------loading the data done!------------")    

##2) set the ligand-receptor interaction database
if (opt_data_base == "original"){
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  
  # use a subset of CellChatDB for cell-cell communication analysis
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
} else if (opt_data_base =="manually curated"){
  # extending the CellChat database - NOTE: the extended database can be created step-wise with the code in the accompanying "mouse" script and "forcing" gene names mouse -> human
  CellChatDB <- CellChatDB.human
  celltalk_cellphone = read.delim("celltalkdb_cellphonedb_harmonized_for_cellchat_HUMAN24082023LFB.txt")
  celltalk_cellphone
  # correcting relevant gene <-> protein names
  ppg = read.table("protein_converted_interactions.txt")
  for(i in 1:nrow(celltalk_cellphone)) {
    
    matching_index <- which(ppg$V1 == celltalk_cellphone$ligand[i])
    
    if(length(matching_index) == 1) {
      celltalk_cellphone$ligand[i] <- ppg$V3[matching_index]
      celltalk_cellphone$pathway_name[i] <- ppg$V3[matching_index]
    }
  }
  
  for(i in 1:nrow(celltalk_cellphone)) {
    
    matching_index <- which(ppg$V1 == celltalk_cellphone$receptor[i])
    
    if(length(matching_index) == 1) {
      celltalk_cellphone$receptor[i] <- ppg$V3[matching_index]
      celltalk_cellphone$interaction_name[i] = paste0(celltalk_cellphone$ligand[i], "_", celltalk_cellphone$receptor[i])
      try({rownames(celltalk_cellphone)[i] = celltalk_cellphone$interaction_name[i]})
    }
  }
  CellChatDB$interaction <- rbind(CellChatDB$interaction, celltalk_cellphone) %>% distinct()
  
  #what's in the data base
  showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  
  # standard CellChat preprocessing steps
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
}
rm(CellChatDB)
rm(CellChatDB.use)
print("----------loading the database done!------------")

##3) preprocessing the data
#To infer the cell state-specific communications, we identify over-expressed ligands or receptors 
#in one cell group and then identify over-expressed ligand-receptor interactions if either ligand 
#or receptor is over-expressed. 
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multisession", workers = 20) # do parallel

if (opt_server==TRUE){
   options(future.globals.maxSize = 3 * 8000 * 1024^2)
} else {
   options(future.globals.maxSize = 8000 * 1024^2)
}

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto Protein-Protein-Interaction (Optional: when running it, USER should set 
#`raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
print("----------preprocessing done!------------")

##4) Inference of cell-cell communication network
#CellChat infers the biologically significant cell-cell communication by assigning each interaction 
#with a probability value and peforming a permutation test. CellChat models the probability of 
#cell-cell communication by integrating gene expression with prior known knowledge of the interactions 
#between signaling ligands, receptors and their cofactors using the law of mass action. 

#When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations 
#tend to send collectively stronger signals than the rare cell populations, CellChat can also consider 
#the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE. 
cellchat <- computeCommunProb(cellchat) # optionally set type="truncatedMean", trim ="0.1" --> default is trimean
saveRDS(cellchat, file = paste0(results_path,"cellchat_",subset_condition,"_before_filtering_DB_",opt_data_base,".rds"))

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
print("----------cell cell communication network inference done!------------")

##5) Infer the cell-cell communication at a signaling pathway level
#The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in 
#the slot 'net' and 'netP', respectively.
cellchat <- computeCommunProbPathway(cellchat)
print("----------cell cell communication network inference at signaling pathway level done!------------")

##6) Calculate the aggregated cell-cell communication network
#by counting the number of links or summarizing the 
#communication probability. USER can also calculate 
#the aggregated network among a subset of cell groups 
#by setting sources.use and targets.use. 
cellchat <- aggregateNet(cellchat)
print("----------calculating aggregated cell cell communication network done!------------")

##7) Extract the inferred cellular communication network as a data frame
#data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
head(df.net)
#inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5. 
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.signalingpws <- subsetCommunication(cellchat,slot.name = "netP" )
#save data frames
saveRDS(df.net, file=paste0(results_path,"df_net",subset_condition,".Rda"))
saveRDS(df.signalingpws, file=paste0(results_path,"df_signalingpws",subset_condition,".Rda"))
print("----------extract + save inferred cellular network done!------------")


##8) save the cellChat object
saveRDS(cellchat, file = paste0(results_path,"cellchat_",subset_condition,"_after_filtering_DB_",opt_data_base,".rds"))
print("----------save cellchat object done!------------")


