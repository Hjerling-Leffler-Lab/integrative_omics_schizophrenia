#author: Lisa Bast, rewritten from loaddata_filtered.R and pagoda_filtering.R (author: José A. Martínez López)
#date: 24.09.2021
#version: 0.0.2
#about: Read matrices from cellranger 6.0.1 (introns included) using Pagoda2
#       perform filtering using Pagoda2, creates Fig. S1A
#       export to loom file with 2 layers (raw and filtered)
#       

library(pagoda2)
library(stringr)

#define paths:
code_path = getwd()
source("utils.R")

setwd("../")
main_path = getwd()
setwd("../")
main_project_path = getwd()
main_data_path <- paste0(main_project_path,"/2_alignment/output/")

libraries <- paste0(replicate(5,"Library_"),seq(1,5))
lib_paths <- paste0(paste0(replicate(5,paste0(main_data_path)),libraries),"/")
n_files_per_library <- c(10,18,18,20,19)
file_id_start <- 1

setwd(main_data_path)
matrix_path<-"/outs/filtered_feature_bc_matrix"

path_output<-paste0(main_project_path,"/3_quality_control/output/")

#list_files_str = ""
#sample_counter <- 1
sample_counter <- 10+18+18+20+file_id_start
sample_names = c()
library_names = c()
sample_names_vec = c()
library_names_vec = c()
list_files <- vector(mode = "list", length = 1)
#list_files <- vector(mode = "list", length = sum(n_files_per_library))
for (l_id in seq(1,length(libraries))){
  for (file_id in seq(file_id_start,n_files_per_library[l_id])){
    sample_names <- append(sample_names,paste0("S",sample_counter))
    library_names <- append(library_names,paste0(libraries[l_id]))
    path_i <- paste0(paste0(paste0(paste0(paste0(paste0(lib_paths[l_id]),"Counts_10x"),str_sub(libraries[l_id], 9)),"_"),file_id),matrix_path)

    #add library and sample name for each cell in current sample 
    new_line <- paste0(paste0(paste0(paste0("S",sample_counter),"='"),path_i),"'")
    list_files[[1]] <- noquote(new_line)

    #read data for current sample:
    counts<-read.10x.matrices(eval(parse(text=list_files)))
    #perform pagoda2 filtering
    sample_str <- paste0(paste0(paste0(paste0("Counts_10x"),str_sub(libraries[l_id], 9)),"_"),file_id)
    tryCatch({ pdf(file = paste0(paste0(paste0(path_results_pagoda_filt,"pagoda_S"),sample_counter),".pdf"),width=10,height=10)
      counts_filtered <- gene.vs.molecule.cell.filter(counts,max.cell.size = 50000,min.cell.size=500)
      dev.off()
    },
    error = function(c) paste0(paste0("error: ",sample_str)," could not be pagoda filtered!"),
    warning = function(c) paste0(paste0("warning: ",sample_str)," could not be pagoda filtered!"),
    message = function(c) paste0(paste0("message: ",sample_str)," could not be pagoda filtered!")
    )
    
    #save as loom file and add sample name and library to metadata
    n_cells_i <- dim(counts)[2]
    sample_names_vec <- rep(paste0("S",sample_counter),n_cells_i)
    library_names_vec <- rep(paste0(libraries[l_id]),n_cells_i)
    tryCatch(save_as_loom_file(counts, counts_filtered, path_output, paste0(paste0("S",sample_counter),"_cellranger_pagoda_filtered_R.loom"), sample_names_vec, library_names_vec),
             error = function(c) paste0(paste0("error: ",sample_str)," could not be saved!"),
             warning = function(c) paste0(paste0("warning: ",sample_str)," could not be saved!"),
             message = function(c) paste0(paste0("message: ",sample_str)," could not be saved!")
    )
    sample_counter <- sample_counter+1
  }
}

