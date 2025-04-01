#about: cluster integrated GO terms with package rrvgo
#       -longread
#       -proteomics
#       -DEGs any cell type
#       -TWAS
#       -any module
#       -any bordeaux module
#       -SCZ risk enriched DEGs
#author:Lisa Bast
#date: 31.05.2024
#version: 0.0.3
# To do: 
# - implement for longread and proteomics
# - run for all the pathways found in any of the analysis to identify parentTerms and export mapping of GO terms to their parent term to use in all kinds of plots later on 


#load libraries
library(rrvgo)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
#library(GOSemSim)

path_code = getwd()
setwd("../")
main_path_project  = getwd()
setwd("../")
main_path  = getwd()

#settings:
opt_input_data_vec = c("DEG_GSA_result_all_CTs","longread_results","TWAS_result","all_GRN_modules_gprofiler_result","all_bordeaux_modules_gprofiler_result","proteomics_result")
opt_input_data_separately = FALSE#TRUE#

#comp_name = 'longleaf'
p_val_cutoff <- 0.05
opt_sim_based_on <- 'size' # "scores"

data_and_results_path <- paste0(main_path,'/output/GSA_Integrated_analysis/')
if (!dir.exists(data_and_results_path)){
  dir.create(data_and_results_path)
} 

setwd(path_code)
source('pathway_clustering_utils.R')

D <- get_all_data(opt_input_data_vec, main_path, p_val_cutoff)
write.csv(D,paste0(data_and_results_path,"pathway_results_all_analyses.csv"))

#extract group, subgroup, geneset, 1-P.bonf.group as score

for (GO_category in c('BP','MF','CC')){
  #filter data set for relevant rows (pathways)
  D_sel <- D %>% filter(subgroup==paste0("GO:",GO_category))
  #browser()
  if (all(colnames(D_sel)!="ID")){
    D_sel$ID <- D_sel %>% {str_extract(.$geneset, "^.{10}")}
  }
  if (opt_input_data_separately==TRUE){
    groups <- unique(D_sel$TestVar)
  } else{
    groups <- "all"
  }
  for (group in groups){
    if (opt_modules_separately ==TRUE){
      D_sel_gr <- D_sel %>% filter(grepl(TestVar, group))
    } else{
      D_sel_gr <- D_sel # nothing filtered, put all modules together
    }
    
    if (dim(D_sel_gr)[1]>2){
      #load_OrgDb("org.Hs.eg.db")
      
      # 1) calculate score similarity matrix
      #browser()
      simMatrix <- calculateSimMatrix(D_sel_gr$ID,orgdb="org.Hs.eg.db",ont=GO_category,method="Rel")
      
      if (!is.null(dim(simMatrix))){
        if (opt_sim_based_on == 'size'){
          reducedTerms <- reduceSimMatrix(simMatrix,
                                          threshold=0.7,
                                          orgdb="org.Hs.eg.db")
        } else if (opt_sim_based_on == 'scores'){
          scores <- setNames(-log10(D_sel_gr$P.bonf.group), D_sel_gr$ID)
          reducedTerms <- reduceSimMatrix(simMatrix,
                                          scores,
                                          threshold=0.7,
                                          orgdb="org.Hs.eg.db")
        }
        #save table of reduced terms:
        filename_RT <- paste0('reduced_terms_rrvgo_',group,'_',GO_category,'_',opt_sim_based_on,'.csv')
        write.csv(reducedTerms,paste0(data_and_results_path,filename_RT))
      }
    }
  }
}

