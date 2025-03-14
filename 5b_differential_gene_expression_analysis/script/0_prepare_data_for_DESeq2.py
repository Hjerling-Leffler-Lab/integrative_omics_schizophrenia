# about: prepare data for DEG (Differentially expressed genes) analysis scRNAseq/ snRNAseq data with DESeq2
# date: Created on Wed Aug  4 13:42:20 2021
# @author: Lisa Bast
# version: 0.0.2

#to do: fix paths and filenames and delete unccessary stuff

import os
import numpy as np
import loompy
import pandas as pd
import sys
import gc

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_main = path_code.replace("5b_differential_gene_expression_analysis\\script","")

opt_pagoda_filtering = True
#opt_celltype_groups = 'scmap'
opt_celltype_groups = 'conos_cluster_based'

opt_aggregate_across = "cell types" #"all" #
opt_df_cells_per_donor_and_CT = 'load' #'create'#
opt_downsample_cells = False#
opt_use_downsampled_loom_files = False
if opt_downsample_cells or opt_use_downsampled_loom_files:
    n_datasets = 10

aggregation_methods = ['sum']#,'mean']#,'median','sum']

add_str = '_cellranger'
n_clusters_range = [15]

if opt_pagoda_filtering:
    add_str_pagoda = '_pagoda'
else:
    add_str_pagoda = ''


for n_clusters in n_clusters_range:
    # read data as anndata object
    path_data,path_results = ut.get_paths(path_main,'DEG_ana_preparation',opt_aggregate_across) 
    
    if opt_aggregate_across=="cell types":
        path_data = "D:/Documents/projects/SCZ_human_cortex/data/CT_clustered_loom_formatted_data_cellranger/15_CTs/"
        
    cl_str = '_'+str(n_clusters)+'_CTs'
    os.chdir(path_data)

    str_not_contains = 'downsampled'
    str_start = 'Samples_cellranger'+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups
    str_end = '.loom'
    loom_files = [ f for f in os.listdir( os.curdir ) if os.path.isfile(f) and f.startswith(str_start) and f.endswith(str_end) and str_not_contains not in f]
    
    CTs_assigned = [f.replace(str_start+"_","").replace(str_end,"") for f in loom_files]
    for opt_aggregation in aggregation_methods:
        #aggregate across cell types
        if opt_aggregate_across == "cell types":
            
            for c_id,ct in enumerate(CTs_assigned):
                #replace '-' in CT string
                ct_str_aggr_data = ct.replace('-','_')  
                ct_loom_filename = 'Samples_cellranger'+add_str_pagoda+ '_TH_and_D_adj_filtered_'+opt_celltype_groups+'_'+ct+'.loom'
                
                ct_csv_matrix_filenames = 'Agg_counts_'+add_str_pagoda+'_TH_and_D_adj_filtered_'+ct_str_aggr_data+'_'+opt_aggregation+'.csv'
                ct_csv_group_info_filenames = 'G_info_'+add_str_pagoda+'_TH_and_D_adj_filtered_'+ct_str_aggr_data+'_'+opt_aggregation+'.csv'
                ut.create_aggregated_data(path_results,path_main,path_data,ct_loom_filename,opt_aggregation, ct_csv_matrix_filenames, ct_csv_group_info_filenames,opt_downsample_cells,False)
        
        #aggregate across all cells 
        if opt_aggregate_across == "all":
            path_data = "D:/Documents/projects/SCZ_human_cortex/data/filtered_loom_formatted_data_cellranger/"
            loom_filename_all_CTs = "Samples_conos_cellranger_pagoda_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"
            ct_csv_matrix_filename_all = 'Agg_counts_'+add_str_pagoda+'_TH_and_D_adj_filtered_all_'+opt_aggregation+'.csv'
            ct_csv_group_info_filename_all = 'G_info_'+add_str_pagoda+'_TH_and_D_adj_filtered_all_'+opt_aggregation+'.csv'
            ut.create_aggregated_data(path_results,path_main,path_data,loom_filename_all_CTs,opt_aggregation, ct_csv_matrix_filename_all, ct_csv_group_info_filename_all,opt_downsample_cells,False)

