# author: Lisa Bast
# date: 2023-11-17,  15:18:21
# version: 0.0.1
# about: clustering modules based on genes they comprise
#        Fig. 5F

import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn import cluster
import seaborn as sns

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project_main = path_code.replace("script","")
path_main = path_code.replace("5c_gene_co_expression_network_analysis\\script","")
path_data = path_project_main + "/4_data_integration_and_cell_type_annotation/output/"
path_module_genes = path_project_main + "/output/cluster_name_15CTs/scz/cor_bicor/"
path_DEG_result = path_main + "/5b_differential_gene_expression_analysis/output/DEGs/16_CTs/design_with_6_RUVs/"
DEG_result_filename = "df_DESeq2_all_genes_all_CTs.csv"
path_genelists = path_main + "/5b_differential_gene_expression_analysis/data/"
gene_matrix_file = path_genelists + "geneMatrix_v2.tsv"

#settings:
n_cl = 15
opt_test_DEM = "LR"
#opt_filter_for_bordeaux_modules = True
alpha_val_DIUGs = 0.1
alpha_val_DEPs = 0.05
#paths
bordeaux_modules = ["Astrocytes_module_yellow",
            "Excitatory_Layer_5-6_IT_neurons_I_module_green",
            "Oligodendrocytes_module_red",
            "Inhibitory_LAMP5_neurons_module_red",
            "Inhibitory_VIP_neurons_module_pink",
            "Excitatory_Layer_2-3_IT_neurons_I_module_red",
            "Excitatory_Layer_2-3_IT_neurons_II_module_black", 
            "Inhibitory_PVALB_neurons_module_magenta"] 

CT_list = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
        #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",

opt_binary_plot=True # changed later for proteomics

sign_modules_names = []
for ct in CT_list:
    DF_module_genes = pd.read_csv(path_module_genes+"module_gene_mapping_info_"+ct+".csv")
    #filter for sign modules:
    print(ct)
    module_names = ut.get_scz_relevant_module_names(ct,path_module_genes,opt_test_DEM,0.05)
    print(module_names)
    ct_mod_names = [ct+'_module_'+m for m in module_names]
    sign_modules_names = sign_modules_names + ct_mod_names
    DF_module_genes_filt = DF_module_genes[DF_module_genes['module_name'].isin(ct_mod_names)]
    DF_module_genes_filt["celltype"] = ct
    if ct==CT_list[0]:
        DF = DF_module_genes_filt.copy()
    else:
        DF = pd.concat([DF,DF_module_genes_filt])
#export sign_module names
np.savetxt(path_module_genes+"sign_module_names.csv",sign_modules_names,delimiter=",",fmt="% s")    
if len(DF)==0:
    next

modes=["count"]
#create matrix of CT_modules vs. pathways with 0s and 1s
for mode in modes:
    M_df, opt_binary_plot = ut.get_matrix_of_modules_vs_pathways(DF,"no_filtering",mode)
    if np.shape(M_df)[0]>0:
        ut.plot_clustermap_modules_vs_genes(M_df,path_module_genes,opt_test_DEM,"no_filtering",opt_binary_plot)
        M_df=M_df.transpose()
        #M = M_df.to_numpy()
        if opt_binary_plot==True: # for counts
            ut.plot_number_genes_per_module(M_df,n_cl,path_module_genes,path_data,opt_test_DEM,"no_filtering") 
            #calculate pairwise correlations of modules:
            DF_C,DF_P = ut.get_module_correlations_and_percentage_overlapping_genes(M_df)

            ##use code understand-gsa-gene-overlap.txt in GSA analyis folder
            module_list_ordered = ut.plot_clustermaps(DF_C,DF_P,M_df,path_module_genes,opt_test_DEM,"no_filtering",n_cl,path_data,gene_matrix_file,opt_binary_plot)
            #export list with order of modules
            df = pd.DataFrame(module_list_ordered, columns = ["module"])
            df.to_csv(path_module_genes+"modules_ordered_"+opt_test_DEM+".csv")

