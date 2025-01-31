# author: Lisa Bast
# date: 2024-06-26,  14:50:41
# version: 0.0.1
# about: put together pathway results of modules in one common excel sheet

import pandas as pd
import os
import sys

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()

main_path= path_code.replace("5c_gene_co_expression_network_analysis\\script","")
path_data = main_path + "/4_data_integration_and_cell_type_annotation/output/"
module_results_path = main_path + "5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/cor_bicor/variable/" 

opt_test = "LR" 
CT_list_modules = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
                #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",

with pd.ExcelWriter(module_results_path+'T8_pathway_results_modules.xlsx', engine='openpyxl') as writer:
    for ct in CT_list_modules:
        folder = module_results_path + ct+"/norm_basic/"
        df = pd.read_csv(folder + "g_profiler_results_"+ct+".csv")
        #columns to keep:
        c_keep = [c for c in df.columns if c not in ["p_values","significant","query_sizes","intersection_sizes","effective_domain_size","source_order"]]
        #save in excel sheets:
        df[c_keep].to_excel(writer, sheet_name=ct)


