# author: Lisa Bast
# date: 2024-04-17,  10:13:00
# version: 0.0.1
# about: a) plot connectivity distribution overall per module
#        b) highlight mitochrondrial DNA genes, RNA binding protein coding genes, lncRNA, other protein coding genes per network module
#        Fig. S20D, S21C
      

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

#paths:
path_project_main = path_code.replace("script","")
path_main = path_code.replace("5c_gene_co_expression_network_analysis\\script","")
path_data = path_project_main + "/4_data_integration_and_cell_type_annotation/output/"
path_module_genes = path_project_main + "/output/cluster_name_15CTs/scz/cor_bicor/"
path_DEG_result = path_main + "/5b_differential_gene_expression_analysis/output/DEGs/16_CTs/design_with_6_RUVs/"
DEG_result_filename = "df_DESeq2_all_genes_all_CTs.csv"
path_genelists = path_main + "/5b_differential_gene_expression_analysis/data/"
gene_matrix_file = path_genelists + "geneMatrix_v2.tsv"
path_results_connectivity = path_module_genes+"/connectivity/"
path_mito_results = path_project_main + "/output/Mitochondria/"
mito_gene_file = "all_nuclear_and_mitochondrial_genes_of_bordeaux_modules.csv"
path_results_proteomics = path_main+"/5b_differential_gene_expression_analysis/data/proteomics/" 
hgnc_file = "hgnc_2022-04-24_small.tsv"
filename_proteomics = "diagnosis_vs_layers.xlsx" #differential_abundance_analysis.xlsx"

#settings:
opt_species = 'human'
opt_velo = False
opt_pagoda_filtering = True
n_cl = 15
opt_with_MT_genes = True
opt_test = "LR"
opt_module_selection = "bordeaux"#"significant_scz" # # "midnight",   "all"# 
opt_metacell_settings = "50_13_250" # "30_6_120"

#gene_matrix_file = path_gene_lists + "geneMatrix_v2.tsv"
gene_ensgid_mapping_file = "Gene_ensgid_mapping.csv"
#path_mBATcombo_results = path_module_genes + "bordeaux_modules/mBATcombo/"

module_CT_list = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
        #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",

modules = ut.get_module_selection(opt_module_selection,module_CT_list,path_module_genes,opt_test)

#load hub genes for module selection:
mapping = pd.read_csv(path_genelists+gene_ensgid_mapping_file)
DF_hub = ut.get_hub_gene_DF(module_CT_list, path_module_genes,modules,mapping) #returns for selected modules only
#get short module names
DF_hub = ut.get_short_CT_names(DF_hub,"celltype")
DF_hub["module_name_short"] = DF_hub["celltype_short"]+'_'+DF_hub["module"]
DF_hub["module_name_short"] =DF_hub["module_name_short"].str.replace("_"," ")


#add proteomics status
#load proteomics data:
DF_p = ut.get_DEP_results(path_results_proteomics,filename_proteomics,1 ,"","", 1 ,path_results_proteomics+hgnc_file)
DF_p["proteomics"]=True
#only significant DEPs
#DF_p = DF_p[DF_p["qvalue_DEP"]<0.05]
DF = pd.merge(left = DF_hub, right=DF_p, left_on="gene_name", right_on="Gene",how="left")

#which hubgenes are more or less abundant proteins in schizophrenia?
#plot connectivity distribution per module, colored by DEP log2FC
ut.plot_connectivity_distribution(DF,opt_module_selection,path_results_connectivity,"proteomics")


#add gene type:
GM = pd.read_csv(path_genelists+"geneMatrix_v2.tsv",sep="\t")
DF_hub = pd.merge(left = GM[["ensgid","gene_type"]],right = DF_hub,on="ensgid",how="right")
DF_hub.loc[DF_hub["gene_type"]=="protein_coding","gene_type"] = "other protein coding gene"

#add RBP status:
#get list of RNA binding protein coding genes
df_RBPs = pd.read_csv(path_genelists+"RBPs_homo_sapiens.txt",sep="\t",header=None)#2960 genes, source: http://eurbpdb.gzsys.org.cn/download.html
df_RBPs.rename(columns={0:"ensgid",1:"gene_name"},inplace=True)
df_RBPs["RBP coding gene"]=True
DF_hub = pd.merge(left=df_RBPs[["ensgid","RBP coding gene"]],right=DF_hub,on="ensgid",how="right")
DF_hub.loc[DF_hub["RBP coding gene"]==True,"gene_type"]="RBP coding gene"
DF_hub.drop(columns="RBP coding gene",inplace=True)

#add nuclear mito gene status:
df_mito = pd.read_csv(path_mito_results+mito_gene_file).transpose() # contains MT genes, too
df_mito["nuclear mito gene"] = True
df_mito = df_mito.reset_index()
df_mito.rename(columns={"index":"gene_name"},inplace=True)
df_mito["gene_name"] = df_mito["gene_name"].str.strip("'")
DF_hub = pd.merge(left=df_mito,right=DF_hub, on="gene_name", how="right")
DF_hub.loc[DF_hub["nuclear mito gene"]==True,"gene_type"]="nuclear mito gene"
DF_hub.drop(columns="nuclear mito gene",inplace=True)

#mitochondrial DNA gene:
DF_hub.loc[DF_hub["gene_name"].str.startswith("MT-"),"gene_type"]="MT gene"
#define other:
DF_hub.loc[~DF_hub["gene_type"].isin(["MT gene","other protein coding gene","nuclear mito gene","RBP coding gene"]),"gene_type"]="other"

#add DEG status
DF_DEGs = pd.read_csv(path_DEG_result+DEG_result_filename,delimiter=',')
DF_hub = pd.merge(left=DF_DEGs[["Gene_short","celltype","log2FoldChange","padj"]],right=DF_hub,right_on=["gene_name","celltype"],left_on = ["Gene_short","celltype"],how="right")
DF_hub["regulation status"] = "unchanged"
DF_hub.loc[np.logical_and(DF_hub["padj"]<0.3,DF_hub["log2FoldChange"]<0),"regulation status"] = "downregulated"
DF_hub.loc[np.logical_and(DF_hub["padj"]<0.3,DF_hub["log2FoldChange"]>0),"regulation status"] = "upregulated"
DF_hub.loc[DF_hub["padj"].isna(),"regulation status"] = "not_tested"

#plot connectivity distribution per module, colored by cell type
ut.plot_connectivity_distribution(DF_hub,opt_module_selection,path_results_connectivity, "DEGs")

