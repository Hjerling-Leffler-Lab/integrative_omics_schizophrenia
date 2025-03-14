# create cell type specific/ cluster specific loom files
#after integration of conos result
# Created on Tue Mar  9 14:25:35 2021
# @author: Lisa Bast
# version: 0.1.0

#import anndata
#import scanpy as sc
import os
import numpy as np
import loompy
import sys
#import pandas as pd

opt_celltype_groups ='conos_cluster_based'
#opt_celltype_groups = 'scmap'
opt_pagoda_filtering = True
opt_D_version = ''#'v1'

add_str_align = '_cellranger'

if opt_pagoda_filtering:
    add_str_pagoda = '_pagoda'
else:
    add_str_pagoda=''
    
if opt_celltype_groups == 'scmap':
    n_clusters_range = [76,51]
elif opt_celltype_groups == 'conos_cluster_based':
    n_clusters_range = [15,37]

os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("4_data_integration_and_cell_type_annotation\\script","")


for n_clusters in n_clusters_range:
    # read data as anndata object
    file_path, _, new_file_path = ut.get_paths(path_project,"Cell type-specific files")  
    
    os.chdir(file_path)
    #aD_i = anndata.read_loom('S1_velocyto_TH_and_D_adj_filtered_and_CT_annotated.loom', X_name='matrix') 
    #aD_i.var_names_make_unique()
    
    if opt_celltype_groups == 'scmap':
        ct_col_of_interest = 'CT_ann_ABM_MCA_scmap_cell2cluster_' + str(n_clusters) + 'CTs'
    elif opt_celltype_groups == 'conos_cluster_based':
        ct_col_of_interest = 'cluster_name_'+str(n_clusters)+'CTs'
    
    source_file = 'Samples_'+opt_celltype_groups[0:5]+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated_and_CT_clustered.loom'
    D=loompy.connect(source_file)
    #with loompy.connect(source_file) as D:
    
        #ct_col_of_interest = 'cluster_name'
    CTs_assigned = ut.get_CT_assigned(D.ca,opt_celltype_groups,ct_col_of_interest)
    #CTs_assigned = np.unique(D.ca[ct_col_of_interest])
    #CTs_assigned = ['Inh L1-6 PVALB SCUBE3 or L3-6 PVALB MFI2','Astro L1-6 FGFR3 ETNPPL or L1 FGFR3 MT1G or FOS']
    for ct in CTs_assigned:
        print(ct)
        genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten() # all genes
        cells_kept = np.zeros((1, np.shape(D)[1]), dtype=bool).flatten() # all cells set to False
        id_cells_kept = np.where(D.ca[ct_col_of_interest]==ct)
        for i in id_cells_kept:
            cells_kept[i]=1 # only cells of this celltype
        new_filename = 'Samples'+add_str_align+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups+'_'+ct.replace(" ","_")+'.loom'
        ut.create_loom_file(D, genes_kept, cells_kept,new_file_path,new_filename)
