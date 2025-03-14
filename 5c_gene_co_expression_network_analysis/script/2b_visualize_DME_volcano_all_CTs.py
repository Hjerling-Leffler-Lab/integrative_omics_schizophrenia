# author: Lisa Bast
# date: 2023-11-30,  13:40:29
# version: 0.0.1
# about: DME result for a specific test all CTs together
#        Fig. 5B

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

path_main = path_code.replace("5c_gene_co_expression_network_analysis\\script","")
path_main_project = path_code.replace("script","")

opt_pagoda_filtering = True
n_cl = 15
opt_test = "LR"#"wilcox" (conservative),  "bimod" (not accurate), "LR", "t", "negbinom", "poisson"
opt_metacell_settings = "50_13_250"

#paths
path_data = path_main_project + "/data/"
path_code = path_main_project + "/script/"
path_results = path_main_project + "/output/"
path_module_genes = path_results + "/cluster_name_15CTs/scz/cor_bicor/" + opt_metacell_settings

if opt_metacell_settings == "50_13_250":
    CT_list = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                    "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                    "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                    "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
            #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
elif opt_metacell_settings == "30_6_120":
    CT_list = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                    "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                    "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                    "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I",
                    "Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells"]
            #no modules as cell type populations are too small: "Endothelial_and_mural_cells",,

bool_first = True
for ct in CT_list:
    filename = path_module_genes+ct+"/norm_basic/"+"DMEs_"+opt_test+".csv"
    if os.path.exists(filename):
        DME_results = pd.read_csv(filename)
        DME_results["statistical test"] = opt_test
        if bool_first:
            R_DME = DME_results
            bool_first = False
        else: 
            R_DME = pd.concat([R_DME,DME_results],axis=0)

R_DME["-log10(adj p-val)"] = np.log10(R_DME["p_val_adj"])*(-1)

#integrate ct color in df:
#color_df = fg.get_CT_colors_as_df(15,path_data)
#R_DME = pd.merge(left=R_DME,right=color_df, left_on="level", right_on="celltype",how="left")
#non-sign grey
#R_DME[R_DME["p_val_adj"]>0.05]["color"] = "grey"

ut.plot_volcano_all_modules(R_DME,n_cl,path_data,path_module_genes,opt_test,True)


