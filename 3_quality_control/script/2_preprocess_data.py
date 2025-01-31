# about: preprocesing pipeline snRNAseq data, threshold filtering, Fig. S1 B, C, D
# date: Created on Tue Apr 21 16:13:41 2020
# @author: Lisa Bast
# version: 0.3.1

## specify project path:
## specify project path:
import os
import sys
import gc

os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("3_quality_control\\script","")

## specify settings:
opt_create_and_merge_raw_data_loom_files = False#True# if TRUE loom files per sample are merged in one large loom file containing all samples
opt_pagoda_filtering = True
if opt_pagoda_filtering:
    opt_merge_pagoda_filtered_loom_files = False
    opt_use_pagoda_filtered_data = True# 
else:
    opt_merge_pagoda_filtered_loom_files = False
    opt_use_pagoda_filtered_data = False

opt_use_TH_filtered_data = True #if TRUE the already threshold filtered data is loaded 
opt_use_TH_and_D_filtered_data = True #if TRUE the already threshold  and doublet filtered data is loaded 
opt_remove_doublets = True
if opt_remove_doublets:
    opt_manual_doublet_TH_adjustment = True
opt_D_version = ''#'v2'
opt_stats_raw_data = 'load' #'calculate' #
opt_stats_filt_data = 'calculate' #'load'#
opt_integrate_filtering_status_in_raw_files = False #True #
opt_plot_raw_data = True#False#
opt_plot_filtered_data = True#False
opt_plot_filtering_stats = True#
opt_save = True
if opt_plot_filtered_data:
    opt_violin_plots = False
    opt_spearman_corr_with_QC_metrics = True

path_data_raw, path_code, path_results_raw = ut.get_paths(path_project,opt_pagoda_filtering,'raw',0)
path_data_filt,_,path_results_filt = ut.get_paths(path_project,opt_pagoda_filtering,'preprocessed',0)
path_summary_metrics = path_project + "/3_quality_control/output/"
## specify libraries and samples to drop:
library_nrs = range(1,6) #[1] #
samples_to_drop = ['10x2_19'] #2_19: experimentally messed up (Fatima);
file_str1 = 'TH_filtering'
if opt_remove_doublets:
    if opt_manual_doublet_TH_adjustment:
        file_str2 = 'TH_and_D_adj_'+opt_D_version+'filtering'
    else:
        file_str2 = 'TH_and_D_filtering'
if opt_pagoda_filtering:
    file_str1 = 'pagoda_'+file_str1
    file_str2 = 'pagoda_'+file_str2
    
## (1) get Threshold values:
    #if you want to adapt them change self_fun_get.py
    #make sure to run code again with opt_use_TH_filtered_data = False
#TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_of_cells_expr_a_gene, TH_min_number_genes, TH_min_number_counts, TH_min_number_of_s_counts_per_gene, TH_min_number_of_us_counts_per_gene, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_max_ribosomal_fraction, TH_min_cells_per_sample = fg.get_TH_values(opt_loom_velo)
TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_genes, TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_max_ribosomal_fraction, TH_min_cells_per_sample = ut.get_TH_values()
## (2) get raw data
D_raw = ut.get_raw_data(library_nrs,samples_to_drop,opt_create_and_merge_raw_data_loom_files,path_project,path_code,path_data_raw)
sample_IDs,n_cells,disease_status,sex_status = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D_raw)
 
## (3) calculate some statistics per sample
# - CPG: cells/ counts per gene
# - MCG: mean counts per gene
# - min_MCG: minimum counts per gene
# - GPC: genes per cell
# - RPC: reads per cell
# - MitoRPC: mitochondrial reads per cell
# - PMito: percentage of reads belonging to mitochondrial mRNA per cell
# - RiboRPC: ribosomal reads per cell
# - PRibo: percentage of reads belonging to ribosomal mRNA per cell
# - MinRPC: minimum reads per cell
# - RLD: distribution of read length (1 value in distribution per gene)
MCG, CPG, min_CPG, GPC, CPC, sCPC, usCPC, MCPC, PMito, RCPC, PRibo, mean_counts_per_cell, std_counts_per_cell, RLD, med_gpc, df_stats = ut.get_stats_per_sample(D_raw,path_results_raw,path_code,opt_stats_raw_data,'',opt_save)

D_raw.close()

## (4) plot raw data
if opt_plot_raw_data==True:
    opt_doublets_removed=False
    ut.plot_data([], [], TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_min_number_genes, TH_min_cells_per_sample, path_project, path_summary_metrics, False, False, False, opt_stats_raw_data, opt_save,opt_D_version, opt_pagoda_filtering)
        
## (5) filter raw data based on thresholds and remove doublets (optional)
if opt_pagoda_filtering == True:
    # identify which cells got filtered and store information in raw data loom files as metadata variable
    D_fil, D_fil_CTRL, D_fil_SCZ, df_pagoda_filtering_result = ut.get_pagoda_filtered_data(sample_IDs, disease_status, path_code, path_data_filt, path_results_filt, opt_create_and_merge_raw_data_loom_files, opt_use_TH_filtered_data, opt_use_TH_and_D_filtered_data, opt_plot_filtered_data, TH_min_cells_per_sample,opt_D_version,opt_use_pagoda_filtered_data)
    #update sample_IDs,n_cells,disease_status,sex_status:
    sample_IDs,n_cells,disease_status,sex_status = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D_fil)
if opt_remove_doublets == True: 
    D_fil, D_fil_CTRL, D_fil_SCZ, df_TH_and_D_filtering_result, PD, PD_updated = ut.get_TH_and_D_filtered_data(sample_IDs, disease_status, path_code, path_data_raw, path_data_filt, path_results_filt, opt_pagoda_filtering, opt_create_and_merge_raw_data_loom_files, opt_use_TH_filtered_data, opt_use_TH_and_D_filtered_data, opt_manual_doublet_TH_adjustment, opt_plot_filtered_data, TH_min_cells_per_sample,opt_D_version)
    opt_doublets_removed = True
    opt_TH_filtered = True
    if opt_plot_filtering_stats:
        ut.plot_percentage_cells_discarded(file_str2,path_results_filt,True)
        ut.plot_percentage_cells_discarded(file_str2,path_results_filt,False)
else:
    D_fil, D_fil_CTRL, D_fil_SCZ, df_TH_filtering_result = ut.get_TH_filtered_data(sample_IDs, disease_status, path_code, path_data_filt, path_results_filt, opt_create_and_merge_raw_data_loom_files,opt_use_TH_filtered_data,opt_pagoda_filtering)
    opt_TH_filtered = True
    
if opt_plot_filtering_stats:
    ut.plot_percentage_cells_discarded(file_str1,path_results_filt,True)
    ut.plot_percentage_cells_discarded(file_str1,path_results_filt,False)

sample_IDs_CTRL_fil,n_cells_CTRL_fil,disease_status_CTRL_fil,sex_status_CTRL_fil = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D_fil_CTRL)
sample_IDs_SCZ_fil,n_cells_SCZ_fil,disease_status_SCZ_fil,sex_status_SCZ_fil = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D_fil_SCZ)
sample_IDs_fil,n_cells_fil,disease_status_fil,sex_status_fil = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D_fil)
# ## remove doublets and save ThresHold and Doublets filtered loom files
## create loom file with all filtered CTRL samples
if (opt_plot_filtering_stats and opt_TH_filtered and opt_doublets_removed):
    if opt_integrate_filtering_status_in_raw_files:
        loom_file_name = ut.integrate_filtering_status_in_metadata(path_data_raw,path_data_filt, opt_D_version,opt_pagoda_filtering)
    else:
        loom_file_name = 'Samples_cellranger.loom'
    ut.plot_UMAP_removed_cells(path_results_filt,path_data_raw,loom_file_name)

# if opt_stats_filt_data =='load': 
filtering_str = ut.get_filtering_str(opt_TH_filtered,opt_doublets_removed,opt_manual_doublet_TH_adjustment,opt_D_version, opt_pagoda_filtering)
MCG_f, CPG_f, min_CPG_f, GPC_f, CPC_f, sCPC_f, usCPC_f, MCPC_f, PMito_f, RCPC_f, PRibo_f, mean_counts_per_cell_f, std_counts_per_cell_f, RLD_f, med_gpc_f, df_stats_f= ut.get_stats_per_sample(D_fil,path_results_filt,path_code,opt_stats_filt_data,filtering_str,opt_save)
#data file D_fil contains all filtered samples (filtered according to settings)
if opt_plot_filtered_data==True:
    ut.plot_data(PD, PD_updated, TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_min_number_genes, TH_min_cells_per_sample, path_project, path_summary_metrics, opt_TH_filtered, opt_doublets_removed, opt_manual_doublet_TH_adjustment, opt_stats_filt_data, opt_save, opt_D_version, opt_pagoda_filtering)   
    ut.plot_expression_PCA('_TH_and_D_adj_filtered.loom',file_str1,file_str2,opt_violin_plots,opt_spearman_corr_with_QC_metrics,path_data_filt,path_results_filt)
