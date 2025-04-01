# -*- coding: utf-8 -*-
"""
about: plot volcano plots specifically for genes of certain GO terms across cell types that were significantly enriched for this GO term
       Fig. 3 G-I
date: Created on Thu Sep 29 14:07:17 2022
@author: Lisa Bast
version: 0.0.1
"""

import pandas as pd
import numpy as np
import os
import sys 
import seaborn as sns

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

main_path = path_code.replace("4_data_integration_and_cell_type_annotation\\script","")

path_colors = main_path + "/4_data_integration_and_cell_type_annotation/output/"
path_results_mito = main_path + "/5b_differential_gene_expression_analysis/output/"
path_genes_matrix = path_results_mito + "/Gene_lists_GO_DB/"
path_DEG_results = path_results_mito + "/DEGs/16_CTs/design_with_6_RUVs/"

#read_DEGs:
DEGs = pd.read_csv(path_DEG_results+"T5_DEGs_per_cell_type.csv")
DEGs[["Gene","ENSGID"]] = DEGs["Gene"].str.rsplit(pat=r'_',n=1,expand=True)

#load genes for respective mito GO terms:
GO_terms = ["mitochondrial outer membrane","mitochondrial inner membrane","neuron differentiation","immune response"]
for term in GO_terms:
    df_ensgids = pd.read_csv(path_genes_matrix + "GO_DB/ensgid_GO_CC_"+term+".tsv",sep="\t")
    df_ensgids[term] = True
    #print(term+": "+str(len(df_ensgids)))
    #filter for those genes:
    DEGs_filtered = pd.merge(left = DEGs, right = df_ensgids, on="ENSGID",how="inner")
    #make gene index
    DEGs_filtered.set_index("Gene",inplace=True)
    if "mito" in term:
        ut.print_stats_volcano(DEGs_filtered,term)
        #all CTs together:
        ut.plot_volcano(DEGs_filtered,term,path_results_mito,0.3,True,annotations="small_pvals")
        #all neurons and astrocytes:
        DEGs_filtered_neuro_astro = DEGs_filtered[np.logical_or(np.logical_or(DEGs_filtered["celltype"].str.startswith("Exc"),DEGs_filtered["celltype"].str.startswith("Inh")),DEGs_filtered["celltype"].str.startswith("Astro"))]
        ut.plot_volcano(DEGs_filtered_neuro_astro,term+"_astro_neuro",path_results_mito,0.3,True,annotations="small_pvals")
        for ct in DEGs_filtered["celltype"].unique().tolist():
            DEGs_filtered_ct_i = DEGs_filtered[DEGs_filtered["celltype"]==ct]
            #plot_volucano
            ut.plot_volcano(DEGs_filtered_ct_i,term+"_"+ct,path_results_mito,0.3,True,annotations="small_pvals")
    if term == "immune response":
        DEGs_filtered_micro = DEGs_filtered[DEGs_filtered["celltype"].str.startswith("Micro")]
        ut.print_stats_volcano(DEGs_filtered_micro,term)
        ut.plot_volcano(DEGs_filtered_micro,term+"_microglia",path_results_mito,0.3,True,annotations="small_pvals")
    if term=="neuron differentiation":
        ut.print_stats_volcano(DEGs_filtered,term)
        #all CTs together:
        ut.plot_volcano(DEGs_filtered,term,path_results_mito,0.3,True,annotations="small_pvals")


