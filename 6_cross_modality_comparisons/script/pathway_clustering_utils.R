get_settings_for_input_data<- function(opt_input_data){
  if (opt_input_data == "Patient_clustering_genes_GSA_result"){
    n_clusters_range = c('SCZ_cases')#,'whole_dataset','controls')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('BrainCortexExpressed/','genesDetectedInSnRNAseq/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c('all')
  } else if (opt_input_data =="Patient_clustering_genes_GSEA_result"){
    n_clusters_range = c('controls','SCZ_cases','whole_dataset')
    min_overlap_range <- c(0)
    background_folders <- c('none/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c('all')
  } else if (opt_input_data == "DEG_GSA_result"){ 
    n_clusters_range <- c(15)#,37,3)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('BrainCortexExpressed/','BrainCortexExpressed_protein_coding_only/','genesDetectedInSnRNAseq/','genesDetectedInSnRNAseq_protein_coding_only/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c("Excitatory_Layer_5_6_IT_neurons_I","Excitatory_Layer_5_6_IT_neurons_II",
                   "Excitatory_Layer_5_6_CT_and_NP_neurons","Excitatory_Layer_2_3_IT_neurons_I",
                   "Excitatory_Layer_2_3_IT_neurons_II", "Excitatory_Layer_3_4_IT_neurons","Excitatory_Layer_3_6_IT_neurons",
                   "Astrocytes","Endothelial_and_mural_cells","Microglial_cells","Oligodendrocyte_progenitor_cells","Oligodendrocytes",
                   "Inhibitory_LAMP5_neurons","Inhibitory_PVALB_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons")                     
    #old definition:                                        
    #CT_groups <- c('all','Astrocytes', "Microglial_cells","Oligodendrocyte_progenitor_cells",'Neurons')
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){ 
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('norm_basic/')
    min_overlap_range <- c(0)
    background_folders <- c('none/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- TRUE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
    n_clusters_range <- c(15)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('BrainCortexExpressed/')#,'BrainCortexExpressed_protein_coding_only/','genesDetectedInSnRNAseq/','genesDetectedInSnRNAseq_protein_coding_only/')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if(opt_input_data == "all_GRN_modules_gprofiler_result"|| opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    n_clusters_range <- c(15)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('BrainCortexExpressed/')#,'BrainCortexExpressed_protein_coding_only/','genesDetectedInSnRNAseq/','genesDetectedInSnRNAseq_protein_coding_only/')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "all_GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('GSA/')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- TRUE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if (opt_input_data == "GRN_modules_GSA_result"){
    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    sub_folders <- c('GSA/')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    opt_modules_separately <- FALSE
    CT_groups <- c('all') ## placeholder: this is later used for the different modules
  } else if(opt_input_data =="longread_results"){
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "TWAS_result"){
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    background_folders <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "proteomics_result"){
    background_folders <- c('BrainCortexExpressed/')
    n_clusters_range <- c('dummy')
    sub_folders <- c('')
    min_overlap_range <- c('')
    opt_sim_based_on <- 'size'
    CT_groups <- c('all')
  } else if (opt_input_data == "DEG_GSA_result_all_CTs"){
    n_clusters_range <- c(15)#,37,3)
    sub_folders <- c('alpha3/')
    min_overlap_range <- c(3)#c(3,5,10)
    background_folders <- c('BrainCortexExpressed/')
    opt_sim_based_on <- 'size' # 'scores'
    #'all' for plotting pathways significant in any cell type together
    CT_groups <- c("Excitatory_Layer_5_6_IT_neurons_I","Excitatory_Layer_5_6_IT_neurons_II",
                   "Excitatory_Layer_5_6_CT_and_NP_neurons","Excitatory_Layer_2_3_IT_neurons_I",
                   "Excitatory_Layer_2_3_IT_neurons_II", "Excitatory_Layer_3_4_IT_neurons","Excitatory_Layer_3_6_IT_neurons",
                   "Astrocytes","Endothelial_and_mural_cells","Microglial_cells","Oligodendrocyte_progenitor_cells","Oligodendrocytes",
                   "Inhibitory_LAMP5_neurons","Inhibitory_PVALB_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons")
  }
  if ((str_detect(opt_input_data,"GRN")) && (opt_input_data != "all_GRN_modules_gprofiler_result") && (opt_input_data != "all_bordeaux_modules_gprofiler_result")){
    n_clusters_range <- get_cell_type_clusters_GRN(opt_metacell_settings)
  }
  return(list(n_clusters_range,sub_folders,min_overlap_range,background_folders,opt_sim_based_on,CT_groups))
}

get_GSA_results_as_DF <- function(data_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data){
  #to do: track cell type
  setwd(data_path)
  D <- c()
  if (opt_up_and_down_sep==FALSE){
    bool_first = TRUE
    if (opt_input_data == "Patient_clustering_genes_GSA_result"){
      file_list <- Sys.glob("*_up_or_down_all.csv")
    } else if (opt_input_data == "Patient_clustering_genes_GSEA_result"){
      file_list <- Sys.glob("result_*.csv")
    } else {
      file_list <- Sys.glob("*__all.csv")
    }
    # opt_input_data=="GRN_modules_GSA_result" | opt_input_data=="all_GRN_modules_GSA_result")
    for (file in file_list){
      d <- read.csv2(file) 
      d <- d %>% select(c('group', 'subgroup', 'geneset','P.bonf.group','TestVar')) %>% filter(P.bonf.group<=p_val_cutoff)
      if (bool_first==FALSE){
        D <- rbind(D,d)
      } else{
        D <- d
        bool_first <- FALSE
      }
    }
  } else{
    D$up <- c()
    D$down <- c()
    bool_first = TRUE
    reg_strings <- c('up','down')
    for (reg_str in reg_strings){
      file_list <- Sys.glob(paste0("*",reg_str,"_all.csv"))
      if (opt_input_data == "Patient_clustering_genes_GSA_result"){
        file_list_or <- Sys.glob(paste0("*_or_*"))
        file_list <- file_list[!(file_list %in% file_list_or)]
      } else if (opt_input_data == "Patient_clustering_genes_GSEA_result"){
        print("set opt_up_and_down_sep to FALSE!")
        D = c() 
      }
      for (file in file_list){
        d <- read.csv2(file) 
        d <- d %>% select(c('group', 'subgroup', 'geneset','P.bonf.group','TestVar')) %>% filter(P.bonf.group<=p_val_cutoff)
        if (bool_first==FALSE){
          if (reg_str == 'up'){
            D$up <- rbind(D$up,d)
          } else{
            D$down <- rbind(D$down,d)
          }
        } else{
          if (reg_str == 'up'){
            D$up <- d
          } else{
            D$down <- d
          }
          bool_first <- FALSE
        }
      }
    }
  }
  return(D)
}

get_gprofiler_results_as_DF <- function(data_path, p_val_cutoff){
  file <- Sys.glob(paste0(data_path,"g_profiler_results_*.csv"))
  d <- read.csv(file) 
  module_columns <- grep('_module_',colnames(d),value=TRUE)
  module_columns_short <- module_columns[endsWith(module_columns,'_grey')==FALSE]
  module_columns_super_short <- as.character(lapply(strsplit(module_columns_short, split="_"),tail, n=1))
  d <- d %>% select(c('source','term_id', 'term_name', module_columns_short))
  colnames(d) <- c('source','term_id', 'term_name', module_columns_super_short)
  #filter modules for SCZ sign? 
  scz_sign<- read.csv(paste0(data_path,"DMEs_LR.csv"))
  sign_modules <- scz_sign$module[scz_sign$p_val_adj<=0.05]
  d <- d %>% select(c('source','term_id', 'term_name', sign_modules))
  D <- pivot_longer(d, cols= sign_modules)
  #rename columns
  colnames(D)[colnames(D)=="name"] <- "Module"
  colnames(D)[colnames(D)=="value"] <- "pvalue"
  colnames(D)[colnames(D)=="term_id"] <- "ID"
  colnames(D)[colnames(D)=="source"] <- "subgroup"
  #filter for 
  D <- D %>% filter(pvalue<=p_val_cutoff)
  return(D)
}

get_longread_results <- function(data_path, p_val_cutoff){
  filename <- paste0(data_path,"go_hypergeoPathway_filtE_fdr10.tsv")
  d <- read.table(file=filename, sep = '\t', header = TRUE)
  d <- d %>% select(c('subgroup', 'geneset','Phyper',"BH")) 
  #d_GO_BP <- dplyr::filter(d, grepl("GO:BP",subgroup))
  #d_GO_CC <- dplyr::filter(d, grepl("GO:CC",subgroup))
  #d_GO_MF <- dplyr::filter(d, grepl("GO:MF",subgroup))
  #d_GO_BP_sign <- d_GO_BP %>% filter(Phyper/dim(d_GO_BP)[1] <= p_val_cutoff)
  d$P.bonf.group <- d$BH
  d$ID <- sub(" .*", "", d$geneset)
  d$TestVar <- "longread"
  d <- d[ , -which(names(d) %in% c("Phyper","BH","group"))]
  return(d)
}

get_TWAS_results <- function(data_and_results_path, p_val_cutoff){
  filename <- paste0(data_and_results_path, "df_15CT_TWAS0.3_GSA0.05_GO_group3class.tsv")
  d <- read.table(file=filename, sep = '\t', header = TRUE)
  d <- d %>% select(c('subgroup', 'geneset','P.fdr.group')) 
  d$P.bonf.group <- d$P.fdr.group
  d$ID <- sub(" .*", "", d$geneset)
  d$TestVar <- "TWAS"
  d <- d[ , -which(names(d) %in% c("P.fdr.group"))]
  return(d)
}

get_proteomics_results <- function(data_and_results_path, p_val_cutoff){
  filename <- paste0(data_and_results_path, "proteomics_both_sign.csv")
  d <- read.table(file=filename, sep = ';', header = TRUE)
  d$ID <- sub(" .*", "", d$geneset)
  d$P.bonf.group <- d$P.fdr.group
  d <- d %>% select(c('subgroup', 'geneset','ID','P.bonf.group')) 
  d$TestVar <- "proteomics"
  return(d)
}

get_data <- function(opt_input_data,main_path, data_and_results_path,p_val_cutoff, cf, sf, opt_up_and_down_sep,opt_metacell_settings){
  #"DEG_GSA_result_all_CTs"
  ##"longread_results"
  #"TWAS_result"
  ##"all_GRN_modules_gprofiler_result"
  #"all_bordeaux_modules_gprofiler_result"
  #"proteomics_result"
  
  if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    
    if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
      D_d <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
    }
    sf <- 'norm_basic/'
    n_c_r <- get_cell_type_clusters_GRN(opt_metacell_settings)
    
    bool_first <- TRUE
    for (n_c in n_c_r){
      cf <- paste0(n_c,'/')
      data_and_results_path_GRN <- paste0(data_and_results_path,cf,sf)
      D_m_tibble_tmp <- get_gprofiler_results_as_DF(data_and_results_path_GRN, p_val_cutoff)
      D_m_tibble_tmp <- D_m_tibble_tmp %>% mutate_at(.vars = "Module", ~ paste0(n_c,"_",.x))
      D_m_tibble_tmp <- D_m_tibble_tmp %>% unite("geneset", ID:term_name, sep=" ", na.rm=TRUE)
      if (bool_first==TRUE){
        D_m_tibble <- D_m_tibble_tmp
        bool_first <- FALSE
      } else{
        D_m_tibble <- rbind(D_m_tibble,D_m_tibble_tmp)
      }
    }
    #harmonize dataframes:
    D_m <- as.data.frame(D_m_tibble)
    if (nrow(D_m)){
      colnames(D_m)[colnames(D_m)=="Module"] <- "TestVar"
      D <- D_m
    }
    if (opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result"){
      if (nrow(D_d$up)>0){
        colnames(D_d$up)[colnames(D_d$up)=="P.bonf.group"] <- "pvalue"
        D_d$up <- subset(D_d$up,select= -c(group))
        if (nrow(D)>0){
          D <- rbind(D,D_d$up)
        } else{
          D <- D_d$up
        }
      }
      if (nrow(D_d$down)>0){
        colnames(D_d$down)[colnames(D_d$down)=="P.bonf.group"] <- "pvalue"
        D_d$down <- subset(D_d$down,select= -c(group))
        if (nrow(D)>0){
          D <- rbind(D,D_d$down)
        } else{
          D <- D_d$down
        }
      }
    } 
    if (opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "all_bordeaux_modules_gprofiler_result"){
      
      if (opt_input_data=="all_bordeaux_modules_gprofiler_result"){
        #filter for bordeaux modules:
        D <- D[which(D$TestVar %in% c("Oligodendrocytes_red","Excitatory_Layer_5-6_IT_neurons_I_green",
                                     "Inhibitory_LAMP5_neurons_red","Excitatory_Layer_2-3_IT_neurons_II_black", 
                                     "Astrocytes_yellow","Inhibitory_PVALB_neurons_magenta",
                                     "Excitatory_Layer_2-3_IT_neurons_I_red","Inhibitory_VIP_neurons_pink")), ]
        D$TestVar <- "bordeaux_modules"
      } else{
        D$TestVar <- "modules"
      }
      #D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
      D$P.bonf.group <- D$pvalue
      D <- D[ ,-which(names(D) %in% c("pvalue"))]
      D$ID <- sub(" .*", "", D$geneset)
    }
    
  } else if (opt_input_data == "GRN_modules_gprofiler_result"){ # implementation for gprofiler output
    D <- get_gprofiler_results_as_DF(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "longread_results"){
    D <- get_longread_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "TWAS_result"){
    D <- get_TWAS_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "proteomics_result"){
    D <- get_proteomics_results(data_and_results_path, p_val_cutoff)
  } else if (opt_input_data == "GRN_modules_GSA_result"){
    D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
  } else{ # implementation for GSA output
    D <- get_GSA_results_as_DF(data_and_results_path, p_val_cutoff, opt_up_and_down_sep,opt_input_data)
    if (opt_input_data=="DEG_GSA_result_all_CTs"){
      if (opt_up_and_down_sep==TRUE){
        D$up$TestVar <- "DEGs_up"
        D$down$TestVar <- "DEGs_down"
        D <- rbind(D$up,D$down)
      } else{
        D$TestVar <- "DEGs"
      }
      D <- D[ ,-which(names(D) %in% c("group"))]
      #add GO ID column:
      D$ID <- sub(" .*", "", D$geneset)
    }
  }
  return(D)
}

get_all_data <- function(opt_input_data_vec, main_path,p_val_cutoff){
  bool_first <- TRUE
  for (opt_input_data in opt_input_data_vec){
    #print(opt_input_data)
    opt <- get_settings_for_input_data(opt_input_data)
    n_clusters_range <- opt[[1]]
    sub_folders <- opt[[2]]
    min_overlap_range <- opt[[3]]
    background_folders <- opt[[4]]
    opt_sim_based_on <- opt[[5]]
    CT_groups <- opt[[6]]
    #"DEG_GSA_result_all_CTs"
    #"longread_results"
    #"TWAS_result"
    #"all_GRN_modules_gprofiler_result"
    #"all_bordeaux_modules_gprofiler_result"
    #"proteomics_result"
    if (opt_input_data=="GRN_modules_gprofiler_result"){
      opt_up_and_down_sep = FALSE #
    } else if (opt_input_data=="Patient_clustering_genes_GSA_result" | opt_input_data=="GRN_modules_GSA_result" | opt_input_data=="all_GRN_modules_GSA_result" | opt_input_data =="all_bordeaux_modules_gprofiler_result" ){
      opt_up_and_down_sep = FALSE #
    } else{
      opt_up_and_down_sep = TRUE#FALSE #
    }
    for (n_clusters in n_clusters_range){
      if (opt_input_data %in% c("Patient_clustering_genes_GSA_result","Patient_clustering_genes_GSEA_result")){
        cluster_folder <- paste0(n_clusters,'/no_normalization/batch_correction_for_Library/')#norm_var_p_vals_01_var_scale_power_1/') #norm_var_p_vals_00001/')#
        #if (n_clusters=='controls'){
        sub_folders <- c('factor_1/','factor_2/','factor_3/')
        #} else{
        #  sub_folders <- c('factor_1/','factor_2/','factor_3/','factor_4/','factor_5/','factor_6/')
        #}
      } else if (opt_input_data == "DEG_GSA_result" || opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data =="all_bordeaux_modules_gprofiler_result" || opt_input_data == "DEG_GSA_result_all_CTs"){  
        if (n_clusters==3){
          cluster_folder <- paste0(n_clusters,'_classes/')
        } else{
          cluster_folder <- paste0(n_clusters,'_CTs/')
        }
      } else if (opt_input_data=="GRN_modules_gprofiler_result"){
        cluster_folder <- paste0(n_clusters,'/')
      } else{
        cluster_folder <- ""
      }
      for (sub_folder in sub_folders){
        for (min_overlap in min_overlap_range){
          for (background_folder in background_folders){
            
            #load significant pathways:
            data_and_results_path <- get_results_path(main_path,opt_input_data,cluster_folder,min_overlap,background_folder,sub_folder,opt_metacell_settings)
            D <- get_data(opt_input_data,main_path, data_and_results_path,p_val_cutoff, cluster_folder, sub_folder, opt_up_and_down_sep,opt_metacell_settings)
            #make sure order of columns is harmonized:
            D <- D[ ,c("subgroup","geneset","ID","P.bonf.group","TestVar")]
            if (bool_first==TRUE){
              D_all <- D
              bool_first <- FALSE
            } else{
              D_all <- rbind(D_all,D)
            }
          }
        }
      }
    }
  }
  return(D_all)
}

get_results_path <- function(main_path,opt_input_data,cluster_folder,min_overlap,background_folder,sub_folder,opt_metacell_settings){
  if (opt_input_data == "DEG_GSA_result"|| opt_input_data == "DEG_GSA_and_GRN_modules_gprofiler_result" || opt_input_data == "DEG_GSA_result_all_CTs"){
    data_and_results_path <- paste0(main_path,'results/GSA_analysis/cellranger/v3/',cluster_folder,sub_folder,'min_overlap_',as.character(min_overlap),'/6_RUV_and_0_MARS/',background_folder)
  } else if (opt_input_data == "all_GRN_modules_gprofiler_result" || opt_input_data == "GRN_modules_gprofiler_result"){
    data_and_results_path <- paste0(main_path,'/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/')
  } else if (opt_input_data == "all_bordeaux_modules_gprofiler_result"){
    data_and_results_path <- paste0(main_path,'/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/bordeaux_modules/')
  } else if(opt_input_data =="longread_results"){
    data_and_results_path <- paste0(main_path,'/5b_differential_gene_expression_analysis/data/longreadRNAseq/')
  } else if (opt_input_data == "TWAS_result"){
    #to do: specify!
    data_and_results_path <- paste0(main_path,"...")
  } else if (opt_input_data == "proteomics_result"){
    data_and_results_path <- paste0(main_path,"/5b_differential_gene_expression_analysis/data/proteomics")
  } else if (opt_input_data== "GRN_modules_GSA_result" || opt_input_data == "all_GRN_modules_GSA_result"){
    data_and_results_path <- paste0(main_path,'/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/GSA_result/')
  } 
  #print(data_and_results_path)
  return(data_and_results_path)
}


plot_pathway_clustering <- function(simMatrix, reducedTerms, results_path, add_str, GO_category, reg_str, ct_group,opt_up_and_down_sep,opt_input_data){
  if (opt_input_data != "GRN_modules_gprofiler_result"){
    if (opt_up_and_down_sep==TRUE){
      subfolder = 'up_and_down_sep/'
    }else{
      subfolder = 'up_and_down_together/'
    }
  } else{
    subfolder = ''
  }
  path <- paste0(results_path,"figures/rrvgo/",subfolder)
  #make sure path exists:
  dir.create(path, recursive=TRUE, showWarnings = FALSE)
  
  graphic_filename_1 <- paste0(path,"heatmap_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_1,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(heatmapPlot(simMatrix,
                    reducedTerms,
                    annotateParent=TRUE,
                    annotationLabel="parentTerm",
                    fontsize=6))
  dev.off()
  
  if (dim(reducedTerms)[1]>4){
    graphic_filename_2 <- paste0(path,"scatter_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
    pdf(file = graphic_filename_2,   # The directory you want to save the file in
        width = 6, # The width of the plot in inches
        height = 6) # The height of the plot in inches
    print(scatterPlot(simMatrix, reducedTerms, addLabel=TRUE))
    dev.off()
  }

  
  # module_list <- treemapPlot(reducedTerms,size="score")
  # #change some colors
  # if (ct_group=="all" & reg_str=='down' & GO_category=='BP'){
  #   module_list <- hand_curate_colors(module_list)
  #   browser()
  # }
  graphic_filename_3 <- paste0(path,"treemap_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_3,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(module_list <- treemapPlot(reducedTerms,size="score"))
  dev.off()
  
  graphic_filename_4 <- paste0(path,"wordcloud_",add_str,'_',ct_group,'_',GO_category,'_',reg_str,".pdf")
  pdf(file = graphic_filename_4,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(wordcloudPlot(reducedTerms, min.freq=1, colors="black"))
  dev.off()
  
  return(module_list)
}


get_cell_type_clusters_GRN <- function(opt_metacell_settings){
  if (opt_metacell_settings == "50_13_250" ){
    n_c_r = c("Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
              "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
              "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
              "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I")  
  } else if (opt_metacell_settings == "30_6_120"){
    n_c_r = c("Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
              "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
              "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
              "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I",
              "Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells")  
  }
  return(n_c_r)
}


# hand_curate_colors <- function(module_list){
#   term_to_adapt <- c("mitochondrial electron transport, NADH to ubiquinone", "mitochondrion organization", "ATP metabolic process","cell migration","positive regulation of endothelial cell proliferation")
#   take_color_from_term <- c("ATP biosynthetic process","ATP biosynthetic process","ATP biosynthetic process","neuron differentiation","neuron differentiation")
#   for (i in seq(1,length(term_to_adapt))){
#     if (term_to_adapt[i] %in% module_list$tm$parentTerm){
#       new_value <- unique(module_list$tm$vColor[module_list$tm$parentTerm==take_color_from_term[i]])
#       if (length(new_value)>1){
#         new_value <- new_value[1]
#       }
#       module_list$tm$vColor[module_list$tm$parentTerm==term_to_adapt[i]] = rep(new_value,times=sum(module_list$tm$parentTerm==term_to_adapt[i]))
#     }
#   }
#   return(module_list)
# }
