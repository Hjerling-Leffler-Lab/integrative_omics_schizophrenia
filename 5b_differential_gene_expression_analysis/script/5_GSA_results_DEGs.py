# -*- coding: utf-8 -*-
"""
about: visualize pathways found based on DEGs (number of genes in pathways and p-values of pathways, comparison with proteomics result)
       Fig. 3B,D,E, S9, S10, S11, S12A
Date: Wed Mar 23 17:00:32 2022
@author: Lisa Bast
version: 0.0.4
"""

import os
import sys 
import pandas as pd
import numpy as np
import xlrd

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_subfolder = "5b_differential_gene_expression_analysis"
path_project = path_code.replace(path_subfolder+"\\script","")
path_info = {"path_project":path_project}
path_info["path_colors"] = path_project + "4_data_integration_and_cell_type_annotation/output"
path_info["subfolders"] = ['alpha3'] #['alpha05','alpha1','alpha3','alpha5']
path_info["min_overlap_folders"] = ['min_overlap_3']#,'min_overlap_5','min_overlap_10']#['min_overlap_1', 'min_overlap_5','min_overlap_10']
opt_results_path = 'pathway_visualization'
_, path_info["path_results"] = ut.get_paths(path_info["path_project"],opt_results_path)  

options = {"opt_server" : False,
            "opt_laptop" : False, #True#
            "opt_species" : 'human',
            "opt_pagoda_filtering" : True,
            "opt_save" : True,
            "opt_GSA_version" : 'v3', #'v2' #'v1'
            "opt_celltype_groups" : 'conos_cluster_based',
            "opt_proteomics" : False,#True,#
            "opt_proteomics_up_and_down_sep" : False,
            "opt_longread" : False,
            "opt_up_and_down_sep" : True,#False,# if DEGS up and down were separately tested
            "opt_plot_n_pathways_heatmap" : False,#True,#
            "opt_plot_pathway_statistics" : False,#True,#
            "opt_plot_PCA_pathways" :  False,#True, 
            "opt_plot_upsets" : False, #True, #
            "opt_aggregate_geneset_pvals_by_averaging" : False,#True#
            "pval_col_to_extract" : 'P.fdr.group',#'P.fdr.all'
            "index_cols" : ['group','subgroup','geneset'],
            "opt_load_dataframe": False,#True,#
            "opt_create_DF_sorted_filtered" : True,#False,#
            "opt_filter_stats_for_GOs": False #True set to False to plot also SYNGO etc
            }

n_cluster=[16]#, 37,3]

#DESeq2_aggregation_methods = ['sum']
DESeq2_design_modes = ['6_RUVs']

n_TH_genes_per_pathway = [10]#[1,5,10] # number of genes a pathway must have in at least 1 cell type in order to be plotted in case alpha is 0.3
background_data_transcriptomics = ['cortex']#,'cortex_pc']#,'SnRNAseqDet','SnRNAseqDet_pc']#'cortex','cortex_pc',
if options["opt_proteomics"]:
    background_data_proteomics = ['cortex']#,'cortex_pc']#,'genesMappingToProteinsDetectedInProteomics','genesMappingToProteinsDetectedInProteomics']
else:
    background_data_proteomics = []

if options["opt_GSA_version"]=='v1':
    variables_groups = ['go','pathway']
    variable_n_genes = 'xVar'
    variable_gene_names = 'nVar_ensgid'
else:
    variables_groups = ['subgroup','group']
    variable_n_genes = 'overlap.TestVar.geneset'#'xVarX'
    variable_gene_names = 'genes.in.geneset.tested'
    
options["opt_plot_rrvgo"] = True #False #if true, pathways grouped into categories by rrvgo will be plotted instead of ungrouped GO terms
if options["opt_plot_rrvgo"]:
    options["opt_plot_rrvgo_modes"] = ['average_for_each_parent','significant_sorted_by_parent']

if options["opt_aggregate_geneset_pvals_by_averaging"]:
   options[" p_val_cutoff"] = 0.2
else:
   options["p_val_cutoff"] = 0.05#2

for n_cl in n_cluster:
    for sf in path_info["subfolders"]:
        alpha_val_DEGs = int(sf.split('alpha')[1])/10
        for ssf in path_info["min_overlap_folders"]:
            for i,dm in enumerate(DESeq2_design_modes):
                if n_cl == 3:
                    path_info["folder_p"] = path_info["path_results"]+"DEGs/"+options["opt_GSA_version"]+'/'+str(n_cl)+'_classes/'+sf+'/'+ssf+'/'+dm+'/'
                else:
                    path_info["folder_p"] = path_info["path_results"]+"DEGs/"+options["opt_GSA_version"]+'/'+str(n_cl)+'_CTs/'+sf+'/'+ssf+'/'+ dm+'/'
                os.chdir(path_info["folder_p"])
                for j,bg in enumerate(background_data_transcriptomics):
                    folders, results_path_proteomics, bg_str = ut.get_folders(options,path_info,background_data_proteomics,sf,ssf,j,bg)
                    for directory in folders:
                        if directory[-1]!='/':
                            directory = directory+'/'
                        path_info["path_results_subdirectory"] = directory
                        if options["opt_plot_pathway_statistics"] or options["opt_create_DF_sorted_filtered"] or options["opt_plot_n_pathways_heatmap"]==True or options["opt_plot_pathway_statistics"]==True or options["opt_plot_PCA_pathways"]==False or options["opt_plot_upsets"]==True:
                            DF_sorted_filtered = ut.generate_GSA_output_and_statistics(path_info, dm, options, variable_gene_names, n_cl,alpha_val_DEGs,variables_groups,variable_n_genes,background_data_proteomics,directory,j,bg_str,results_path_proteomics)
                        else:
                            DF_sorted_filtered = ut.get_DF_sorted_filtered([],[],[],options["p_val_cutoff"],options["index_cols"],options["opt_create_DF_sorted_filtered"], options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"], background_data_proteomics, j, path_info["path_results_subdirectory"], options["opt_up_and_down_sep"])
                        #get short cell type names as columns
                        DF_sorted_filtered = ut.get_short_CT_names_as_columns(DF_sorted_filtered)
                        ut.plot_triangular_heatmaps(DF_sorted_filtered,path_info, options,n_cl,bg_str)
                        


