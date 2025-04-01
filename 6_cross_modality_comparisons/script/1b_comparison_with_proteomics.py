# author: Lisa Bast
# date: 2023-11-07,  15:08:09
# version: 0.0.1
# about: compare proteomics with snRNA-seq
#           a) pathways in proteomics with module pathways
#           b) module genes coding for DAPs
#        Fig. S21 A,B, S22A

import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()

#paths:
path_project_main = path_code.replace("script","")
path_main = path_code.replace("6_cross_and_inter_modality_comparisons\\script","")
path_DEG_result = path_main + "/5b_differential_gene_expression_analysis/output/DEGs/16_CTs/design_with_6_RUVs/"
DEG_filename = "df_DESeq2_all_genes_all_CTs.csv"
path_genelists = path_main + "/5b_differential_gene_expression_analysis/data/"
gene_matrix_file = path_genelists + "geneMatrix_v2.tsv"
path_filtered_data = path_main + "/4_data_integration_and_cell_type_annotation/output/"
path_filtered_data_per_CT = path_filtered_data + "/CT_specific_files/16_CTs/"
path_results_proteomics_DAPs = path_main+"/5b_differential_gene_expression_analysis/data/proteomics/" 
hgnc_file = "hgnc_2022-04-24_small.tsv"
filename_proteomics_DAPs = "diagnosis_vs_layers.xlsx" #differential_abundance_analysis.xlsx"
path_pw_results_proteomics = "/5b_differential_gene_expression_analysis/output/GSA_analysis/DAPs/"
file_proteomics_pws = path_pw_results_proteomics + "proteomics_both_sign.csv"
path_module_genes = path_main + "/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/"
results_path = path_project_main + "/output/"
#settings:
opt_test_DEM = "LR"
n_cl = 16
alpha_val_DAPs = 0.05

module_CT_list = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
        #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",

bordeaux_modules = ["Astrocytes_module_yellow",
            "Excitatory_Layer_5-6_IT_neurons_I_module_green",
            "Oligodendrocytes_module_red",
            "Inhibitory_LAMP5_neurons_module_red",
            "Inhibitory_VIP_neurons_module_pink",
            "Excitatory_Layer_2-3_IT_neurons_I_module_red",
            "Excitatory_Layer_2-3_IT_neurons_II_module_black", 
            "Inhibitory_PVALB_neurons_module_magenta"] 

#load helper functions:
sys.path.append(path_code)
import utils as ut

##gene level comparison:
opt_binary_plot=True # changed later for proteomics
sign_modules_names = []

for ct in module_CT_list:
    DF_module_genes = pd.read_csv(path_module_genes+"module_gene_mapping_info_"+ct+".csv")
    #filter for sign modules:
    print(ct)
    module_names = ut.get_scz_relevant_module_names(ct,path_module_genes,opt_test_DEM,0.05)
    print(module_names)
    ct_mod_names = [ct+'_module_'+m for m in module_names]
    sign_modules_names = sign_modules_names + ct_mod_names
    DF_module_genes_filt = DF_module_genes[DF_module_genes['module_name'].isin(ct_mod_names)]
    DF_module_genes_filt["celltype"] = ct
    if ct==module_CT_list[0]:
        DF = DF_module_genes_filt.copy()
    else:
        DF = pd.concat([DF,DF_module_genes_filt])
#export sign_module names
np.savetxt(path_module_genes+"sign_module_names.csv",sign_modules_names,delimiter=",",fmt="% s")    
    
_,_,_,DF = ut.get_DEP_DEG_module_dataframe(DF,path_genelists,path_module_genes,path_DEG_result, DEG_filename,path_results_proteomics_DAPs,filename_proteomics_DAPs,alpha_val_DAPs,path_results_proteomics_DAPs,hgnc_file,bordeaux_modules)

opt_filter_genes = "DAP_genes"+str(alpha_val_DAPs)[2:]

modes = ["log2FoldChange_DAP","count"]

#create matrix of CT_modules vs. pathways with 0s and 1s
for mode in modes:
    M_df, opt_binary_plot = ut.get_matrix_of_modules_vs_pathways(DF,opt_filter_genes,mode)
    if np.shape(M_df)[0]>0:
        ut.plot_clustermap_modules_vs_genes(M_df,results_path,opt_test_DEM,opt_filter_genes,opt_binary_plot)
        M_df=M_df.transpose()
        #M = M_df.to_numpy()
        if opt_binary_plot==True: # for counts
            ut.plot_number_genes_per_module(M_df,n_cl,results_path,path_module_genes,opt_test_DEM,opt_filter_genes) 
            #calculate pairwise correlations of modules:
            DF_C,DF_P = ut.get_module_correlations_and_percentage_overlapping_genes(M_df)

            module_list_ordered = ut.plot_clustermaps(DF_C,DF_P,M_df,results_path,opt_test_DEM,opt_filter_genes,n_cl,path_module_genes,gene_matrix_file,opt_binary_plot)
            #export list with order of modules
            if opt_filter_genes=="no_filtering":
                df = pd.DataFrame(module_list_ordered, columns = ["module"])
                df.to_csv(path_module_genes+"modules_ordered_"+opt_test_DEM+".csv")

## pathway level comparison:
#get module pathway results:
bool_first = True
for go_str in ["BP","CC","MF"]:
    file_sn = path_module_genes+"DF_sign_modules_"+go_str+"_all_cell_types_"+opt_test_DEM+".csv"
    if bool_first == True:
        df_sn = pd.read_csv(file_sn)
        bool_first = False
    else:
        df_sn_tmp = pd.read_csv(file_sn)
        df_sn = pd.concat([df_sn,df_sn_tmp])

#get proteomics pathway results:
df_pr = pd.read_csv(file_proteomics_pws,sep=';')
df_pr = df_pr[df_pr["subgroup"].isin(["GO:BP","GO:CC","GO:MF"])].copy()
df_pr_sign = df_pr[df_pr['P.fdr05.group']].copy()


#determine overlap:
df_pr_sign[["go","term"]] = df_pr_sign["geneset"].str.split(" ",n=1,expand=True)
df_overlap = pd.merge(left = df_pr_sign[["go","P.fdr.group"]], right = df_sn, on="go",how="left")
df_overlap.dropna(subset=["Celltype_Module"],inplace=True)

#short module names:
df_overlap = ut.get_short_CT_names(df_overlap,"Celltype_Module")
df_overlap['Celltype_Module_short'] = df_overlap['Celltype_Module_short'].str.replace('_',' ')

for opt_visualize_pvalue in [True,False]:
    if opt_visualize_pvalue==False:
        df_overlap["counts"] = 1
        df_overlap_p = pd.pivot(df_overlap,columns = "Celltype_Module_short", index="term",values="counts")
        df_overlap_p.fillna(0, inplace=True)
    else:
        df_overlap["P.fdr.group"] = df_overlap["P.fdr.group"].str.replace(",",".").astype("float")
        df_overlap_p = pd.pivot(df_overlap,columns = "Celltype_Module_short", index="term",values="P.fdr.group")
        df_overlap_p.fillna(1, inplace=True)
    #make sure all modules are in table
    modules_ordered = pd.read_csv(path_module_genes+"modules_ordered_"+opt_test_DEM+".csv")
    modules_ordered.drop(columns="Unnamed: 0", inplace=True)
    for m in modules_ordered['module'].tolist():
        if m not in df_overlap_p.columns:
            df_overlap_p[m]=0
    df_overlap_p = df_overlap_p[modules_ordered['module'].tolist()]
    if opt_visualize_pvalue==False:
        df_overlap_p.fillna(0, inplace=True)
    else:
        df_overlap_p.fillna(1, inplace=True)
    #plot Figure S21A
    ut.plot_clustermap_module_pathway_comparison(df_overlap_p,True,path_module_genes,"proteomics",opt_test_DEM,opt_visualize_pvalue)
    ut.plot_clustermap_module_pathway_comparison(df_overlap_p,False,path_module_genes,"proteomics",opt_test_DEM,opt_visualize_pvalue)


### log2FC of mitochondrial DNA genes/ proteins:
#get data frame
DEG_MT = ut.get_MT_genes_from_DEG_and_DAP(path_DEG_result,DEG_filename,path_results_proteomics_DAPs,filename_proteomics_DAPs,hgnc_file,path_filtered_data)
#plot Figure S22B
ut.plot_MT_genes_log2FC_per_CT(DEG_MT,results_path)

