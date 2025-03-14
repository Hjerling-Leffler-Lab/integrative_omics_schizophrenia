# author: Lisa Bast
# date: Wed Apr 22 11:22:02 2020
# version: 0.0.1
# about: helper functions required for preparing data and DEG results visualization

import os
import numpy as np
import loompy
import pandas as pd
import seaborn as sns
from sctriangulate.colors import build_custom_continuous_cmap
import anndata
import scanpy as sc
import math
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
import matplotlib.ticker as tkr
import matplotlib.colors as mcolors
from matplotlib.tri import Triangulation
from upsetplot import plot, from_indicators
import gc
from bioinfokit import visuz
from sklearn.decomposition import PCA
import fnmatch
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from itertools import groupby
import re

FS=16
G=['CTRL','SCZ']
G2=['m','f']

def add_gene_names(DF):
    l=DF['Unnamed: 0'].tolist()
    DF['genes'] = [i.split('_', 1)[0] for i in l]
    return DF

def clear_ws():
#Clears all the variables from the workspace of the spyder application.
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue
     
        del globals()[var]

def get_folders(options,path_info,background_data_proteomics,sf,ssf,j,bg):
    
    if options["opt_proteomics"]:
        if bg == background_data_proteomics[j]:
            bg_str = '(background '+bg+')'
        else:
            bg_str =  '(background '+bg+'/ '+background_data_proteomics[j]+')'
    else:
        bg_str = '(background '+bg+')'
    folders = [path_info["folder_p"]+bg+'/']
    if options["opt_proteomics"]==True:
        results_path_proteomics = path_info["path_project"] + 'output/GSA_analysis/proteomics/v3/'+sf+'/'+ssf+'/'+background_data_proteomics[j]+'/'
    else:
        results_path_proteomics = ''
        
    return folders, results_path_proteomics, bg_str

### create functions:
def create_loom_file(D, genes_kept, cells_kept, path,filename):
    #row_attrs, col_attrs = get_col_and_row_attributes_as_dictionary_for_selection_of_cells_and_genes(D, genes_kept, cells_kept
    row_attrs = { "Gene": D.ra["Gene"][genes_kept],
                      "Accession": D.ra["Accession"][genes_kept]}
    #do this in a more automated way --> failed so far
    #also need: "Conos_cluster_res_5":Di.ca["Conos_cluster_res_5"],
    #"Conos_cluster_res_7":Di.ca["Conos_cluster_res_7"],
    #"Conos_cluster_res_2nd_level"
    
    col_attrs_names = D.ca.keys()
    col_attrs={"CellID": D.ca["CellID"][cells_kept]}
    for col_a in col_attrs_names:
        if col_a != "CellID":
            col_attrs[col_a] = D.ca[col_a][cells_kept]
            

    os.chdir(path)
    M_Sel_rows = np.compress(genes_kept,D[:,:],axis=0)
    M_Sel = np.compress(cells_kept,M_Sel_rows,axis=1)
    loompy.create(filename, M_Sel, row_attrs, col_attrs)
    
def create_aggregated_data(path_results,path_main,path_data,ct_loom_filename,opt_aggregation,  ct_csv_matrix_filenames, ct_csv_group_info_filenames,opt_downsample_cells,opt_mock_bulk):
    current_path = os.getcwd()
    if opt_downsample_cells:
        if opt_mock_bulk:
            subfolder_str = 'Mock_bulk_data/'
        else:
            subfolder_str = 'Downsampled_data/'
    else:
        subfolder_str = ''
    #aggregate data:
    df_m,df_g = get_aggregated_data_frames(ct_loom_filename, opt_aggregation, path_data,path_main)
    #save aggregated data:
    df_m.to_csv(path_results+subfolder_str+ct_csv_matrix_filenames,index=False)
    df_g.to_csv(path_results+subfolder_str+ct_csv_group_info_filenames,index=False)
    os.chdir(current_path)

### get functions:
def get_paths(path_project,opt_results_path,opt_aggregate_across):
    if opt_results_path=='DEG_ana_preparation':
        if opt_aggregate_across=="cell types":
            path_data = path_project + "/4_data_integration_and_cell_type_annotation/output/CT_specific_files/"
            path_results = path_project + '/5b_differential_gene_expression_analysis/output/CT_clustered_aggregated_data/'
        if opt_aggregate_across=="all":
            path_data = path_project + '/4_data_integration_and_cell_type_annotation/output/'
            path_results = path_project + '/5b_differential_gene_expression_analysis/output/aggregated_data/'
    elif opt_results_path.startswith('DEG_visualization'):
        path_results = path_project + '/5b_differential_gene_expression_analysis/output/DEGs/'
        path_data = path_project + '/5b_differential_gene_expression_analysis/data/CT_clustered_aggregated_data/'
    elif opt_results_path.startswith('pathway_visualization'):
        path_results = path_project + '/5b_differential_gene_expression_analysis/output/GSA_analysis/'
        path_data = path_project + '/5b_differential_gene_expression_analysis/output/'
    if not os.path.exists(path_results):
        os.makedirs(path_results)
    return path_data, path_results 

def add_gene_names(DF):
    l=DF['Unnamed: 0'].tolist()
    DF['genes'] = [i.split('_', 1)[0] for i in l]
    return DF

def get_ensemble_genes(df_DESeq2_all_genes_all_CTs,alpha_val,folder,opt_per_CT):
    if opt_per_CT:
        if 'celltype_name' in df_DESeq2_all_genes_all_CTs.columns:
            column_CT = 'celltype_name' 
        else:
            column_CT = 'celltype' 
        CTs = np.unique(df_DESeq2_all_genes_all_CTs[column_CT]).tolist()
        for CT_i in CTs:
            genes_sign_down = df_DESeq2_all_genes_all_CTs[(df_DESeq2_all_genes_all_CTs[column_CT]==CT_i) & (df_DESeq2_all_genes_all_CTs['padj']<=alpha_val) & (df_DESeq2_all_genes_all_CTs['log2FoldChange']<0)]['Gene'].tolist()
            genes_sign_up = df_DESeq2_all_genes_all_CTs[(df_DESeq2_all_genes_all_CTs[column_CT]==CT_i) & (df_DESeq2_all_genes_all_CTs['padj']<=alpha_val) & (df_DESeq2_all_genes_all_CTs['log2FoldChange']>0)]['Gene'].tolist()
            # export list of DEGs's ENSEMBLE IDs:
            genes_sign_up_ensemble = [i.split('_')[1] for i in genes_sign_up]
            genes_sign_down_ensemble = [i.split('_')[1] for i in genes_sign_down]
            genes_sign_up_names = [i.split('_')[0] for i in genes_sign_up]
            genes_sign_down_names = [i.split('_')[0] for i in genes_sign_down]
            DF_up =pd.DataFrame({'genes_ensemble_up':genes_sign_up_ensemble,'genes_names_up':genes_sign_up_names})
            DF_down =pd.DataFrame({'genes_ensemble_down':genes_sign_down_ensemble,'genes_names_down':genes_sign_down_names})
            #create directory if does not exist yet
            path=folder+'genes_ensemble/'
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            DF_up.to_csv('genes_ensemble'+str(alpha_val).split('.')[1]+'_up_'+ CT_i +'.csv', sep=',', columns = ['genes_ensemble_up'], header=False, index=False)
            DF_down.to_csv('genes_ensemble'+str(alpha_val).split('.')[1]+'_down_'+ CT_i +'.csv', sep=',', columns = ['genes_ensemble_down'], header=False, index=False)
            DF_up.to_csv('genes_names'+str(alpha_val).split('.')[1]+'_up_'+ CT_i +'.csv', sep=',', columns = ['genes_names_up'], header=False, index=False)
            DF_down.to_csv('genes_names'+str(alpha_val).split('.')[1]+'_down_'+ CT_i +'.csv', sep=',', columns = ['genes_names_down'], header=False, index=False)
    else:     
        genes_sign_down = df_DESeq2_all_genes_all_CTs[(df_DESeq2_all_genes_all_CTs['padj']<=alpha_val) & (df_DESeq2_all_genes_all_CTs['log2FoldChange']<0)]['Gene'].tolist()
        genes_sign_up = df_DESeq2_all_genes_all_CTs[(df_DESeq2_all_genes_all_CTs['padj']<=alpha_val) & (df_DESeq2_all_genes_all_CTs['log2FoldChange']>0)]['Gene'].tolist()
        genes_intersect_up_down = list(set(genes_sign_down) & set(genes_sign_up))
        # export list of DEGs's ENSEMBLE IDs:
        genes_sign_up_ensemble = [i.split('_')[1] for i in genes_sign_up]
        genes_sign_down_ensemble = [i.split('_')[1] for i in genes_sign_down]
        genes_sign_up_names = [i.split('_')[0] for i in genes_sign_up]
        genes_sign_down_names = [i.split('_')[0] for i in genes_sign_down]
        DF_up =pd.DataFrame({'genes_ensemble_up':genes_sign_up_ensemble,'genes_names_up':genes_sign_up_names})
        DF_down =pd.DataFrame({'genes_ensemble_down':genes_sign_down_ensemble,'genes_names_down':genes_sign_down_names})
        #create directory if does not exist yet
        path=folder+'genes_ensemble/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
            
        os.chdir(path)
        DF_up.to_csv('genes_ensemble'+str(alpha_val).split('.')[1]+'_up.csv', sep=',', columns = ['genes_ensemble_up'], header=False, index=False)
        DF_down.to_csv('genes_ensemble'+str(alpha_val).split('.')[1]+'_down.csv', sep=',', columns = ['genes_ensemble_down'], header=False, index=False)
        DF_up.to_csv('genes_names'+str(alpha_val).split('.')[1]+'_up.csv', sep=',', columns = ['genes_names_up'], header=False, index=False)
        DF_down.to_csv('genes_names'+str(alpha_val).split('.')[1]+'_down.csv', sep=',', columns = ['genes_names_down'], header=False, index=False)
        if len(genes_intersect_up_down)>0:
            DF_intersect =pd.DataFrame({'genes_intersect_up_down':genes_intersect_up_down})
            DF_intersect.to_csv('genes_intersect_up_down'+str(alpha_val).split('.')[1]+'.csv', sep=',', columns = ['genes_intersect_up_down'], header=False, index=False)

def save_DF(DF_sorted_filtered,opt_transcriptomics_up_and_down_sep,opt_proteomics,folder,background_data_proteomics,j,add_str):
    if opt_transcriptomics_up_and_down_sep:
        #make sure same cell types exist in up and down 
        CT_list_up = DF_sorted_filtered.columns[DF_sorted_filtered.columns.str.endswith('_up')]
        CT_list_down = DF_sorted_filtered.columns[DF_sorted_filtered.columns.str.endswith('_down')]
        if len(CT_list_up) < len(CT_list_down):
            #if ct is in down but not in up copy column from _down and replace all values with 1s
            for ct in CT_list_down:
                if ct.split(r'_down')[0]+'_up' not in CT_list_up:
                    DF_sorted_filtered[ct.split(r'_down')[0]+'_up'] = np.nan
                    print(ct.split(r'_down')[0]+'_up'+' added!')
        elif len(CT_list_down) < len(CT_list_up):
            #if ct is in down but not in up copy column from _down and replace all values with 1s
            for ct in CT_list_up:
                if ct.split(r'_up')[0]+'_down' not in CT_list_down:
                    DF_sorted_filtered[ct.split(r'_up')[0]+'_down'] = np.nan
                    print(ct.split(r'_up')[0]+'_down'+' added!')
        
        #export dataframe to be able to load later on for different analysis e.g. mito analysis
        if opt_proteomics:
            DF_sorted_filtered.to_csv(folder+'All_CTs_up_and_down_'+background_data_proteomics[j]+'as_bg_4_proteomics'+add_str+'.csv', header=True, index=True)
        else:
            DF_sorted_filtered.to_csv(folder+'All_CTs_up_and_down'+add_str+'.csv', header=True, index=True)
    else:
        DF_sorted_filtered.to_csv(folder+'All_CTs_together'+add_str+'.csv', header=True, index=True)
    return DF_sorted_filtered

def get_DF_sorted_filtered(DF_up,DF_down,DF_proteomics_all,p_val_cutoff,index_cols, opt_create_DFs,opt_proteomics,opt_proteomics_up_and_down_sep, background_data_proteomics, j, folder,opt_transcriptomics_up_and_down_sep):
    if opt_create_DFs:
        if opt_transcriptomics_up_and_down_sep:
            #put together both dataframes:
            DF_all = pd.merge(DF_down, DF_up,how='outer',on=index_cols)
        else:
            DF_all = DF_up
        
        if (opt_proteomics and (opt_proteomics_up_and_down_sep==False)):
            DF_all = pd.merge(DF_all, DF_proteomics_all,how='outer',on=index_cols)
            
        DF_sorted = DF_all.sort_values(index_cols)
        num_cols = [x for x in DF_sorted.columns if x not in index_cols]
        DF_sorted = save_DF(DF_sorted,opt_transcriptomics_up_and_down_sep,opt_proteomics,folder,background_data_proteomics,j,"_not_filtered")
        DF_sorted_filtered = DF_sorted.loc[DF_sorted[num_cols].min(numeric_only=True,axis=1)<=p_val_cutoff,:]
        DF_sorted_filtered = save_DF(DF_sorted_filtered,opt_transcriptomics_up_and_down_sep,opt_proteomics,folder,background_data_proteomics,j,"")
        
    else: #load data frame instead:
        if opt_transcriptomics_up_and_down_sep:
            if opt_proteomics:
                DF_sorted_filtered = pd.read_csv(folder+'All_CTs_up_and_down_'+background_data_proteomics[j]+'as_bg_4_proteomics.csv')
            else:
                DF_sorted_filtered = pd.read_csv(folder+'All_CTs_up_and_down.csv')
        else:
            DF_sorted_filtered = pd.read_csv(folder+'All_CTs_together.csv')

    return DF_sorted_filtered

def get_short_CT_names_as_columns(DF):
    CTs = DF.columns.tolist()
    CTs_short = [sub.replace("Inhibitory","Inh") for sub in CTs]
    CTs_short = [sub.replace("Excitatory","Exc") for sub in CTs_short]
    #if np.any([" " in c for c in CTs])==False: # underscore instead of space
    CTs_short = [sub.replace("2_3","2-3") for sub in CTs_short]
    CTs_short = [sub.replace("3_4","3-4") for sub in CTs_short]
    CTs_short = [sub.replace("5_6","5-6") for sub in CTs_short]
    CTs_short = [sub.replace("3_6","3-6") for sub in CTs_short]
    CTs_short = [sub.replace("_"," ") for sub in CTs_short]
    CTs_short = [sub.replace("Layer_","L ") for sub in CTs_short]
    CTs_short = [sub.replace("Layer","L") for sub in CTs_short]
    CTs_short = [sub.replace("L ","L") for sub in CTs_short]
    CTs_short = [sub.replace("neurons_","") for sub in CTs_short]
    CTs_short = [sub.replace(" neurons","") for sub in CTs_short]
    CTs_short = [sub.replace("_neurons","") for sub in CTs_short]
    CTs_short = [sub.replace('Oligodendrocyte_progenitor_cells',"OPCs") for sub in CTs_short]
    CTs_short = [sub.replace('Oligodendrocyte progenitor cells',"OPCs") for sub in CTs_short]
    CTs_short = [sub.replace('Microglial_cells',"Microglia") for sub in CTs_short]
    CTs_short = [sub.replace('Microglial cells',"Microglia") for sub in CTs_short]
    CTs_short = [sub.replace('mural_cells',"mural") for sub in CTs_short]
    CTs_short = [sub.replace('mural cells',"mural") for sub in CTs_short]
    CTs_short = [sub.replace("_and_","/") for sub in CTs_short]
    CTs_short = [sub.replace(" and ","/") for sub in CTs_short]
    CTs_short = [sub.replace("  "," ") for sub in CTs_short]
    
    col_rename_dict = {i:j for i,j in zip(CTs,CTs_short)}
    DF.rename(columns=col_rename_dict, inplace=True)
    return DF

def get_aggregated_data_frames(ct_loom_filename, opt_aggregation, path_data, path_main):
    print(path_data)
    os.chdir(path_data)
    with loompy.connect(ct_loom_filename) as D_CT:
        Donor_list = np.unique(D_CT.ca['Donor'])
        bool_first = True
        for donor in Donor_list:
            IDs = np.where(D_CT.ca['Donor']==donor)
            if opt_aggregation == 'sum':
                C_aggr_Si = np.sum(D_CT[:,IDs[0]],axis=1)
            elif opt_aggregation== 'mean':
                C_aggr_Si = np.mean(D_CT[:,IDs[0]],axis=1)
            elif opt_aggregation == 'median':
                C_aggr_Si = np.median(D_CT[:,IDs[0]],axis=1)
            if bool_first:
                df_m = pd.DataFrame(data={'Gene': D_CT.ra['Gene'],'Accession': D_CT.ra['Accession'], donor: C_aggr_Si})
                group_of_donors = [np.unique(D_CT.ca['Disease'][IDs[0]])[0]]
                age_of_donors = [np.unique(D_CT.ca['Age'][IDs[0]])[0]]
                library_of_donors = [np.unique(D_CT.ca['Library'][IDs[0]])[0]]
                pmi_of_donors = [np.unique(D_CT.ca['PMI_h'][IDs[0]])[0]]
                sex_of_donors = [np.unique(D_CT.ca['Sex'][IDs[0]])[0]]
                # p_cell_del_filt = [np.unique(D_CT.ca['p_cell_del_filt'][IDs[0]])[0]]
                # if np.isnan(np.unique(D_CT.ca['p_cell_del_filt'][IDs[0]])[0]):
                #     print('mom mal!')
                # num_reads = [np.unique(D_CT.ca['num_reads'][IDs[0]])[0]]
                # num_umi = [np.unique(D_CT.ca['num_umi'][IDs[0]])[0]]
                # mean_reads_per_umi = [np.unique(D_CT.ca['mean_reads_per_umi'][IDs[0]])[0]]
                bool_first = False
            else:
                df_m[donor] = C_aggr_Si
                group_of_donors.append(np.unique(D_CT.ca['Disease'][IDs[0]])[0])
                age_of_donors.append(np.unique(D_CT.ca['Age'][IDs[0]])[0])
                library_of_donors.append(np.unique(D_CT.ca['Library'][IDs[0]])[0])
                pmi_of_donors.append(np.unique(D_CT.ca['PMI_h'][IDs[0]])[0])
                sex_of_donors.append(np.unique(D_CT.ca['Sex'][IDs[0]])[0])
                # p_cell_del_filt.append(np.unique(D_CT.ca['p_cell_del_filt'][IDs[0]])[0])
                # if np.isnan(np.unique(D_CT.ca['p_cell_del_filt'][IDs[0]])[0]):
                #     print('mom mal!')
                # num_reads.append(np.unique(D_CT.ca['num_reads'][IDs[0]])[0])
                # num_umi.append(np.unique(D_CT.ca['num_umi'][IDs[0]])[0])
                # mean_reads_per_umi.append(np.unique(D_CT.ca['mean_reads_per_umi'][IDs[0]])[0])
                
        df_g = pd.DataFrame(data={'donor_ID_python': Donor_list, 'Group': group_of_donors, 'Age': age_of_donors,'Library':library_of_donors,'PMI':pmi_of_donors,'Sex':sex_of_donors}, index=Donor_list)
        
        os.chdir(path_main+'/3_quality_control/output/')
        qc_metrics = pd.read_excel('T2_quality_metrics_per_sample.xlsx')
        df_g = pd.merge(qc_metrics[['donor_ID_python','p_cells_del_filt','num_reads','num_umi','mean_reads_per_umi']],df_g,on='donor_ID_python')

    return(df_m,df_g)

def get_DEP_results(path_results_proteomics,filename_proteomics,TH_qval_DEPs_padj_DEGs,tissue_layer_str,method_str, alpha_val,hgnc_file):
    if filename_proteomics == 'results_differential_abundance_analysis_proteomics_L4_5_6.xlsx':
        DF_p_ns = pd.read_excel(path_results_proteomics+filename_proteomics,sheet_name = 'Non-signif',engine='openpyxl')
        DF_p_s = pd.read_excel(path_results_proteomics+filename_proteomics,sheet_name = 'All-signif',engine='openpyxl')
        DF_p = pd.concat([DF_p_s,DF_p_ns])
        DF_p_sel = DF_p[['protein_id', 'qvalue', 'foldchange.log2','gene_symbols_or_id']].copy() 
        #rename columns:
        DF_p_sel.rename(columns={'gene_symbols_or_id':'Gene','qvalue':'qvalue_DEP', 'foldchange.log2':'log2FoldChange_DEP'},inplace=True)
    elif filename_proteomics == "diagnosis_vs_layers.xlsx":
        DF_p = pd.read_excel(path_results_proteomics+filename_proteomics,engine='openpyxl')
        DF_p_sel = DF_p[['protein_id', 
                         'qvalue_diagnosis',  
                         'log2fc_diagnosis', 
                         'hgnc_symbol',
                         'hgnc_id']].copy() 
        DF_p_sel.rename(columns={'qvalue_diagnosis':'qvalue_DEP','log2fc_diagnosis':'log2FoldChange_DEP','hgnc_symbol':"Gene"},inplace=True)
    else:
        DF_p = pd.read_excel(path_results_proteomics+filename_proteomics,sheet_name = 'statistics',engine='openpyxl')
        DF_p_sel = DF_p[['protein_id', 
                         'qvalue_'+method_str+'_contrast: control_'+tissue_layer_str+' vs scz_'+tissue_layer_str, 
                         'foldchange.log2_'+method_str+'_contrast: control_'+tissue_layer_str+' vs scz_'+tissue_layer_str, 
                         'gene_symbols_or_id']].copy() 
        #rename columns:
        DF_p_sel.rename(columns={'gene_symbols_or_id':'Gene',
                                'qvalue_'+method_str+'_contrast: control_'+tissue_layer_str+' vs scz_'+tissue_layer_str:'qvalue_DEP', 
                                'foldchange.log2_'+method_str+'_contrast: control_'+tissue_layer_str+' vs scz_'+tissue_layer_str:'log2FoldChange_DEP'},inplace=True)
    DF_p_sel['Color_sign_DEP'] = DF_p_sel['qvalue_DEP'] <= alpha_val
    #gene_symbols_or_id can be separated by ;
    # filter for qval<TH_qval_DEPs_padj_DEGs:
    DF_p_sel = DF_p_sel.loc[DF_p_sel['qvalue_DEP']<=TH_qval_DEPs_padj_DEGs]
    DF_p_sel.loc[(DF_p_sel.Color_sign_DEP==True),'Color_sign_DEP']='dodgerblue'
    DF_p_sel.loc[(DF_p_sel.Color_sign_DEP==False),'Color_sign_DEP']='grey'

    #add ensgid through hgnc:
    DF_hgnc = pd.read_csv(hgnc_file,sep='\t')
    DF_hgnc["HGNC ID"] = "HGNC:"+DF_hgnc["HGNC ID"].astype(str)
    DF_hgnc_sel = DF_hgnc[["HGNC ID", "Ensembl gene ID", "Locus group"]].copy()

    DF_p_sel = pd.merge(left = DF_p_sel, right = DF_hgnc_sel, left_on='hgnc_id', right_on='HGNC ID', how='left')

    return DF_p_sel

def get_color_maps_for_anndata(aD):
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_76CTs_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_51_CTs_colors"] = 'nipy_spectral'
    aD.uns["Disease_colors"] = ['#487BAF','#FA9835']#['steelblue','darkorange']
    aD.uns["Sex_colors"] = ['#CB5C93','#0080FF']#['palevioletred','dodgerblue']
    aD.uns["filtering_status_colors"] = ["#CC0066","#078383","#939393","#2284C6","#000000"]
    return aD

def get_comparison_DEG_pvals_to_genelists(df_DESeq2_all_genes,opt_save,path_results_long,CT_i,alpha_val,cutoff_q_p,bool_first_CT,df_corr_genelists_up,df_corr_genelists_down,path_project,genelist_version):
    
    genelists, genelist_names, qvalues = get_gene_lists_with_qvalues(path_project, genelist_version)

    corr_up = np.empty((len(genelists),1))
    corr_down = np.empty((len(genelists),1))
    corr_up[:] = np.NaN
    corr_down[:] = np.NaN
    n_datapoints_up = np.zeros((len(genelists),1))
    n_datapoints_down = np.zeros((len(genelists),1))
    for i_gl,gl in enumerate(genelists):
        corr_up[i_gl],corr_down[i_gl],n_datapoints_up[i_gl],n_datapoints_down[i_gl] = plot_correlation_DEG_pvalues_with_genelist_p_or_qvalues(df_DESeq2_all_genes,gl,genelist_names[i_gl],qvalues[i_gl],opt_save,path_results_long,CT_i,alpha_val,cutoff_q_p)
        #TO DO: save corr in dataframe and save for each dm and am under path_result_long
        #gene list name vs CT folder name
    if bool_first_CT:
        df_corr_genelists_up = pd.DataFrame(np.hstack([corr_up, n_datapoints_up]),columns =[CT_i+'_corr',CT_i+'_datapoints'] , index=genelist_names)
        df_corr_genelists_down = pd.DataFrame(np.hstack([corr_down, n_datapoints_down]),columns =[CT_i+'_corr',CT_i+'_datapoints'] , index=genelist_names)
        bool_first_CT = False
    else:
        df_corr_genelists_up[CT_i+'_corr'] = corr_up.flatten()
        df_corr_genelists_up[CT_i+'_datapoints'] = n_datapoints_up.flatten()   
        df_corr_genelists_down[CT_i+'_corr'] = corr_down.flatten()
        df_corr_genelists_down[CT_i+'_datapoints'] = n_datapoints_down.flatten()  
    return(df_corr_genelists_up,df_corr_genelists_down,bool_first_CT)

def get_ordered_heatmap_data_genelist_DEG_overlap_and_number_of_genes(DF_red, genelist_name, n_round, opt_par, opt_genes_to_show,opt_q_val_in_gl,mode,opt_add_proteomics):
    DF_red['Gene'] = DF_red.index
    index_col = 'Gene'
    if mode == 'p-value':
        DF_red['padj_range'] = 1-DF_red['padj']
        DF_red.loc[DF_red['log2FoldChange']<0,'padj_range'] =  DF_red.loc[DF_red['log2FoldChange']<0,'padj_range']*(-1)
        var_of_interest = 'padj_range'
    elif mode=='log2FC':
        var_of_interest = 'log2FoldChange'
    if 'celltype_short' in DF_red.columns:
        heatmap_data = pd.pivot_table(DF_red, values = var_of_interest, index=index_col, columns='celltype_short',  dropna=False)
    elif 'celltype_name' in DF_red.columns:
        heatmap_data = pd.pivot_table(DF_red, values = var_of_interest, index=index_col, columns='celltype_name',  dropna=False)
    else:
        heatmap_data = pd.pivot_table(DF_red, values = var_of_interest, index=index_col, columns='celltype', dropna=False)
    if opt_genes_to_show=='top_genes':
        #order columns
        if opt_add_proteomics:
            cols_ordered = heatmap_data.columns[heatmap_data.columns.str.startswith('Prot')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Exc')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Inh')].tolist() + heatmap_data.columns[np.logical_and(np.logical_and(~heatmap_data.columns.str.startswith('Exc'), ~ heatmap_data.columns.str.startswith('Inh')), ~ heatmap_data.columns.str.startswith('Prot'))].tolist() #heatmap_data.columns[1:len(heatmap_data.columns)-1].tolist()+[heatmap_data.columns[0]]
        else:
            cols_ordered = heatmap_data.columns[heatmap_data.columns.str.startswith('Exc')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Inh')].tolist() + heatmap_data.columns[np.logical_and(np.logical_and(~heatmap_data.columns.str.startswith('Exc'), ~ heatmap_data.columns.str.startswith('Inh')),~ heatmap_data.columns.str.startswith('Prot'))].tolist() #heatmap_data.columns[1:len(heatmap_data.columns)-1].tolist()+[heatmap_data.columns[0]]
        heatmap_data = heatmap_data[cols_ordered].copy()
        number_of_genes_plotted = np.shape(heatmap_data)[0]
        number_of_genes_total = 'nan'
    elif opt_genes_to_show=='genes_with_strong_signal':
        number_of_genes_total = np.shape(heatmap_data)[0]
        #remove all rows (genes) that would be completely white anyways (all p-values smaller 0.5, reverse this filtering step as p-values are tranformed to match color bar legend)
        heatmap_data = heatmap_data[np.logical_and(heatmap_data.min(axis=1)>((-1)*(1-opt_par)),heatmap_data.max(axis=1)<(1-opt_par))==False]
        number_of_genes_plotted = np.shape(heatmap_data)[0]
        #order columns
        if opt_add_proteomics:
            cols_ordered = heatmap_data.columns[heatmap_data.columns.str.startswith('Prot')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Exc')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Inh')].tolist() + heatmap_data.columns[np.logical_and(~heatmap_data.columns.str.startswith('Exc'), ~ heatmap_data.columns.str.startswith('Inh'))].tolist() #heatmap_data.columns[1:len(heatmap_data.columns)-1].tolist()+[heatmap_data.columns[0]]
        else:
            cols_ordered = heatmap_data.columns[heatmap_data.columns.str.startswith('Exc')].tolist() + heatmap_data.columns[heatmap_data.columns.str.startswith('Inh')].tolist() + heatmap_data.columns[np.logical_and(~heatmap_data.columns.str.startswith('Exc'), ~ heatmap_data.columns.str.startswith('Inh'))].tolist() #heatmap_data.columns[1:len(heatmap_data.columns)-1].tolist()+[heatmap_data.columns[0]]
        #cols_ordered = heatmap_data.columns[1:len(heatmap_data.columns)-1].tolist()+[heatmap_data.columns[0]]
        heatmap_data = heatmap_data[cols_ordered].copy()

    return (heatmap_data,number_of_genes_total, number_of_genes_plotted)  

def get_percentage_wrongly_assigned_group_labels(path_results):
    #true to compare random labels:
    os.chdir(path_results)
    P=[]
    for i in range(1,101):
        DF = pd.read_csv('group_labels_random_'+str(i)+'.csv')
        DF['group_label_disagreement'] = np.where(DF['group_labels_random']==DF['group_labels_ordered'],False,True)
        #DF.loc[DF['group_labels_random']==DF['group_labels_ordered'], 'group_label_disagreement'] = False
        #DF.loc[DF['group_labels_random']!=DF['group_labels_ordered'], 'group_label_disagreement'] = True
        P = np.append(P,DF['group_label_disagreement'].sum())
    P = P*100/len(P)
    return(P)

def get_percentage_genes_in_list_and_DEG_for_range_of_cutoffs(path_project, DF_DEG, opt_par, qval_cutoff_range):   
    genelists = []
    genelist_names = []
    q_values = []
    #get genelist info
    for genelist_version in ['v1','v2']:
        genelists_i, genelist_names_i, q_values_i = get_gene_lists_with_qvalues(path_project,genelist_version)
        genelists.extend(genelists_i)
        genelist_names.extend(genelist_names_i)
        q_values.extend(q_values_i)
    #for each genelist determine number of genes overlapping with genes investigated DEG analysis and number of genes with signal in DEG analysis
    N_genes_total = np.zeros((len(genelists),len(qval_cutoff_range))) 
    N_genes_plotted = np.zeros((len(genelists),len(qval_cutoff_range))) 
    for i,genelist_i in enumerate(genelists):
        #build dataframe together with DEGs
        if genelist_names[i] in ['asd.wes','scz.wes']:
            string1 = 'q-value '
        else:
            string1 = 'p-value '
        DF = get_genelist_DEG_DF(DF_DEG, genelist_i, genelist_names[i], q_values[i],string1)
        for j,qval_cutoff in enumerate(qval_cutoff_range):
            #get number of genes for currnt genelist 
            DF = get_genelist_DEG_DF(DF_DEG, genelist_i, genelist_names[i], q_values[i],string1)
            #reduce to q-values below qval_cutoff:
            DF_red_all = DF[DF[string1+genelist_names[i]]<qval_cutoff].copy()
            DF_red_all.sort_values(by=string1+genelist_names[i],inplace=True)
            _,idx =np.unique(DF_red_all['Gene'], return_index=True)
            genes_of_interest = DF_red_all['Gene'][np.sort(idx)].tolist()            
            DF_red = DF_red_all[DF_red_all['Gene'].isin(genes_of_interest)]
            if np.shape(DF_red)[0]>0:
                _, N_genes_total[i,j], N_genes_plotted[i,j] = get_ordered_heatmap_data_genelist_DEG_overlap_and_number_of_genes(DF_red, genelist_names[i], 7, qval_cutoff, 'genes_with_strong_signal',True,'p-value',False)
    #calculate percentages:
    P_genes = 100*(N_genes_plotted/N_genes_total)
    return P_genes, genelist_names

def get_gene_lists(path_project,version):
    if version!='synaptic genes':
        os.chdir(path_project+'data/Gene_lists/')

    if version == 'v1' or version == 'v2':
        df_gm = pd.read_csv("geneMatrix_"+version+".tsv",sep='\t')
        #df_gm.head()
    if version =='v1':
        list_autism_wes = df_gm[df_gm['asd.wes']==True]['gene_name']
        list_scz_wes = df_gm[df_gm['scz.wes']==True]['gene_name']
        list_devdelay_wes = df_gm[df_gm['devdelay.wes']==True]['gene_name']
        list_intellectual_disability = df_gm[df_gm['CNVid']==True]['gene_name']
    
        list_adhd2018 = df_gm[df_gm['adhd2018.gene.P']<=0.001]['gene_name']
        list_an2019= df_gm[df_gm['an2019.gene.P']<=0.001]['gene_name']
        list_asd2019= df_gm[df_gm['asd2019.gene.P']<=0.001]['gene_name']
        list_bip2019= df_gm[df_gm['bip2019.gene.P']<=0.001]['gene_name']
        list_mdd2019= df_gm[df_gm['mdd2019.gene.P']<=0.001]['gene_name']
        list_scz2021= df_gm[df_gm['scz2021.gene.P']<=0.001]['gene_name']
        list_iq2018 = df_gm[df_gm['iq2018.gene.P']<=0.001]['gene_name']
       
        genelists = [list_autism_wes,list_scz_wes,list_devdelay_wes,list_intellectual_disability,list_adhd2018,list_an2019,list_asd2019,list_bip2019,list_mdd2019,list_scz2021,list_iq2018]
        genelist_names = ['asd.wes','scz.wes','devdelay.wes','list_CNVid', 'list_adhd2018','list_an2019','list_asd2019','list_bip2019','list_mdd2019','list_scz2021','list_iq2018']
    elif version=='v2':
        list_dd2020_wes = df_gm[df_gm['dd2020.wes.fdr05']==True]['gene_name'] #T/F Dev delay WES, Kaplanis Nature 2020
        list_asd2020_wes = df_gm[df_gm['asd2020.wes.fdr05']==True]['gene_name'] #T/F ASD WES, Satterstrom Cell 2020
        list_asd2021_wes = df_gm[df_gm['asd2021.wes.fdr05']==True]['gene_name'] #T/F ASD WES, Fu 2021 bioRxiv
        list_dd2021_wes = df_gm[df_gm['dd2021.wes.fdr05']==True]['gene_name'] #T/F Dev delay WES, Fu 2021 bioRxiv
        list_ndd2021_wes = df_gm[df_gm['ndd2021.wes.fdr05']==True]['gene_name'] #T/F Neuro dev disorder WES, Fu 2021 bioRxiv
        list_scz2022_wes_bonfSig = df_gm[df_gm['scz2022.wes.bonfSig05']==True]['gene_name'] #T/F SCZ WES, SCHEMA, Nature 2022, bonferroni correction
        list_scz2022_wes_qval = df_gm[df_gm['scz2022.wes.qvalSig05']==True]['gene_name'] #T/F SCZ WES, SCHEMA, Nature 2022, qvalue
        
        list_scz2022InLocus =  df_gm[df_gm['scz2022InLocus']==True]['ensgid']#	T/F	pgc scz2022 gwas, gene in LD-defined locus
        list_scz2022eqtl = df_gm[df_gm['scz2022eqtl']==True]['gene_name']#  T/F		pgc scz2022 gwas, eQTL-gene for eSNP in locus (GTEx cortex)
        list_scz2022PriorityGene = df_gm[df_gm['scz2022PriorityGene']==True]['ensgid']#	T/F	 pgc scz2022 gwas, prioritized gene (fine mapping, eQTL, HiC, rare variant)
        list_scz2022hic = df_gm[df_gm['scz2022hic']==True]['gene_name']#	T/F	 pgc scz2022 gwas, Hi-C anchor in locus and anchor near gene TSS (HCRCI)
        
        #synaptic genes:
        list_from_geneset_tsv = pd.read_csv(path_project+'code/GSA_analysis/gsa-files-v2/genesets.tsv',sep='\t')
        #extract syngo genes
        list_synaptic_genes = list_from_geneset_tsv[list_from_geneset_tsv['group']=='synaptic-gene-ontology']['ensgid'].tolist()
        
        #add CNV genes implicated in scz "list_CNV_del_scz"
        list_CNV_del_scz = df_gm[np.logical_and(df_gm['CNVdel']==True, df_gm['CNVdelPheno'].isin(["adis","a__s","ad_s"]))]['ensgid']
        list_CNV_dup_scz = df_gm[np.logical_and(df_gm['CNVdup']==True, df_gm['CNVdupPheno'].isin(["adis","a__s","ad_s","_d_s"]))]['ensgid']

        genelists = [list_dd2020_wes,list_asd2020_wes,list_asd2021_wes,list_dd2021_wes,list_ndd2021_wes,list_scz2022_wes_bonfSig,list_scz2022_wes_qval,list_scz2022InLocus,list_scz2022PriorityGene, list_scz2022eqtl, list_scz2022hic, list_synaptic_genes, list_CNV_del_scz, list_CNV_dup_scz]
        genelist_names = ['list_dd2020.wes','list_asd2020.wes','list_asd2021.wes','list_dd2021.wes','list_ndd2021.wes','list_scz2022.wes.bonfSig05','list_scz2022.wes.qvalSig05','list_scz2022InLocus','list_scz2022PriorityGene','list_scz2022eqtl','list_scz2022hic', 'list_synaptic_genes','list_CNV_del_scz','list_CNV_dup_scz']
    elif version=='mildas_panel':
       df = pd.read_csv("gene_list_scz_mildas_panel.txt",sep=',',header=None)
       df = df[0].str.replace('[','')
       df = df.str.replace(']','')
       df = df.str.replace('\'','')
       list_from_mildas_panel = df.tolist()
       genelists = [list_from_mildas_panel]
       genelist_names = ['spatial transcriptomics scz panel']

    return genelists, genelist_names

def add_CT_colors_to_DF(DF,path_filtered_loom_data,opt_without_removed_cells,number_clusters):

    colors = get_CT_colors_as_df(number_clusters,path_filtered_loom_data)
    if opt_without_removed_cells == False:
        colors.loc[len(colors)] = {"celltype":"none (removed)","color":DF["celltype_color"][DF["celltype"]=="none (removed)"].unique()[0]}
    if number_clusters==3:
        if "celltype" in DF.columns:
            DF = pd.merge(left=DF, right=colors, on="celltype")
        else:
            colors.rename(columns={"celltype":"cell_type_class","color":"cell_type_class_color"},inplace=True)
            DF = pd.merge(left=DF, right=colors, on="cell_type_class")
    else:
        DF = pd.merge(left=DF, right=colors, on="celltype")
    return DF

def get_CT_colors_as_df(n_clusters,path_filtered_data):
    df = pd.read_excel(path_filtered_data+'/Cell_type_colors.xlsx',engine='openpyxl')
    if n_clusters==3:
        str_ = 'class'
    else:
        str_ = 'type'
    CT_col = 'Cell '+str_+' ('+str(n_clusters)+')'
    Col_col = 'Cell '+str_+' ('+str(n_clusters)+') color'
    _, idx = np.unique(df[CT_col].tolist(),return_index=True)
    DF = df[[CT_col,Col_col]].loc[np.sort(idx)]
    DF=DF.rename(columns={CT_col:'celltype',Col_col:'color'})

    return DF

def get_gene_type_for_deregulated_genes(DF_DEG_results_red, LFC_col, gene_matrix_file,celltype_col,gene_col):
    # to do: distinguish in protein coding and non-protein coding
    #load geneMatrix file
    GM = pd.read_csv(gene_matrix_file,sep='\t')
    #merge genetype info from GM into DF_DEG_result_red
    if gene_col == "Gene_short":
        DF_DEG_results_detailed = pd.merge(left=DF_DEG_results_red,right=GM[["gene_name","gene_type"]],left_on=gene_col,right_on="gene_name", how="left")
    else:
        DF_DEG_results_detailed = pd.merge(left=DF_DEG_results_red,right=GM[["ensgid","gene_type"]],on=gene_col, how="left")
    #decide for each gene if up- or downregulated and add to data frame
    DF_DEG_results_detailed['up-regulated protein coding']=False
    DF_DEG_results_detailed['up-regulated protein coding'][np.logical_and(DF_DEG_results_detailed[LFC_col]>0,DF_DEG_results_detailed["gene_type"]=="protein_coding")]=True
    DF_DEG_results_detailed['down-regulated protein coding']=False
    DF_DEG_results_detailed['down-regulated protein coding'][np.logical_and(DF_DEG_results_detailed[LFC_col]<0,DF_DEG_results_detailed["gene_type"]=="protein_coding")]=True
    DF_DEG_results_detailed['up-regulated lncRNA']=False
    DF_DEG_results_detailed['up-regulated lncRNA'][np.logical_and(DF_DEG_results_detailed[LFC_col]>0,DF_DEG_results_detailed["gene_type"]=="lncRNA")]=True
    DF_DEG_results_detailed['down-regulated lncRNA']=False
    DF_DEG_results_detailed['down-regulated lncRNA'][np.logical_and(DF_DEG_results_detailed[LFC_col]<0,DF_DEG_results_detailed["gene_type"]=="lncRNA")]=True
    DF_DEG_results_detailed['total']=True
    DF_number = DF_DEG_results_detailed.groupby(celltype_col).sum()
    DF_number.drop(columns=[c for c in DF_number.columns if c not in ['up-regulated protein coding','up-regulated lncRNA','down-regulated protein coding','down-regulated lncRNA','total']],inplace=True)

    return DF_number

def get_CT_order(n_clusters,path_filtered_data,opt_short_names):
    if path_filtered_data.endswith('cellranger/'):
        df = pd.read_excel(path_filtered_data+'Cell_type_colors.xlsx',engine='openpyxl')
        if n_clusters==3:
            str_ = 'class'
        else:
            str_ = 'type'
        col_of_interest = 'Cell '+str_+' ('+str(n_clusters)+')'
        if opt_short_names:
            df = get_short_CT_names(df,col_of_interest)
            _, idx = np.unique(df[col_of_interest+"_short"].tolist(),return_index=True)
            CT_order = df[col_of_interest+"_short"][np.sort(idx)].tolist()
        else:
            _, idx = np.unique(df[col_of_interest].tolist(),return_index=True)
            CT_order = df[col_of_interest][np.sort(idx)].tolist()
    else:
        print('This is not implemented yet.')
        CT_order=[]
    return CT_order


def get_number_of_deregulated_genes(DF_DEG_results_red, LFC_col, celltype_col, genes_col, n_clusters, path_filtered_data):

    #decide for each gene if up- or downregulated and add to data frame
    DF_DEG_results_red ['up-regulated']=False
    DF_DEG_results_red ['up-regulated'][DF_DEG_results_red[LFC_col]>0]=True
    DF_DEG_results_red ['down-regulated']=False
    DF_DEG_results_red ['down-regulated'][DF_DEG_results_red[LFC_col]<0]=True
    DF_DEG_results_red ['total']=True
    #for very old pandas version:
    #DF_number = DF_DEG_results_red.set_index([celltype_col,genes_col]).sum(level=celltype_col)[['up-regulated','down-regulated','total']].sort_values(by='total',ascending=False)
    DF_number = DF_DEG_results_red.groupby(celltype_col).sum()[['up-regulated','down-regulated','total']].sort_values(by='total',ascending=False)
    if n_clusters==3:
        new_index=DF_number.index.to_list()
        for i in range(0,len(DF_number)):
            if DF_number.index[i].startswith('NonNeu'):
                new_index[i] = DF_number.index[i]+" cells"
            else:
                new_index[i] = DF_number.index[i]+" neurons"
        DF_number[celltype_col]=new_index
        DF_number.set_index(celltype_col,inplace=True)
    CT_ordered = get_CT_order(n_clusters,path_filtered_data,False)
    #make sure every CT is in DF_number
    if np.shape(DF_number)[0]<len(CT_ordered):
        #CT_missing = [ct.replace(' ','_') for ct in CT_ordered if ct.replace(' ','_') not in DF_number.index.tolist()]
        CT_missing = [ct for ct in CT_ordered if ct not in DF_number.index.tolist()]
        for ct_m in CT_missing:
            DF_number.loc[ct_m]=[0,0,0]
    #order CTs:
   # DF_number = DF_number.reindex(CT_ordered)
    return DF_number

def get_unique_values_in_list_while_preserving_order(input_list):
    _, idx = np.unique(input_list, return_index=True)
    idx_sorted = sorted(idx)
    unique_list_with_preserved_order = [input_list[i] for i in idx_sorted]
    
    return unique_list_with_preserved_order

def get_data_from_df(DF_p,mode,n_clusters, path_filtered_data):
    if mode=='all':
        #if are there columns ending with '_both' or '_all' add them to both and change ending to _up and _down
        valuesNW = DF_p[DF_p.columns[(DF_p.columns.str.endswith('up')|DF_p.columns.str.endswith('_both_x')|DF_p.columns.str.endswith('_all_x'))]]
        valuesSE = DF_p[DF_p.columns[(DF_p.columns.str.endswith('down')|DF_p.columns.str.endswith('_both_y')|DF_p.columns.str.endswith('_all_y'))]]
    elif mode=='proteomics_only':
        if 'proteomics_all' in DF_p.columns:
            valuesNW = DF_p[DF_p.columns[DF_p.columns.str.startswith('proteomics_all')]]
            valuesSE = DF_p[DF_p.columns[DF_p.columns.str.startswith('proteomics_all')]]
        else:
            valuesNW = DF_p[DF_p.columns[DF_p.columns.str.startswith('proteomics_up')]]
            valuesSE = DF_p[DF_p.columns[DF_p.columns.str.startswith('proteomics_down')]]
    elif mode=='transcriptomics_only':
        valuesNW = DF_p[DF_p.columns[(DF_p.columns.str.endswith(' up')&~(DF_p.columns.str.startswith('proteomics')))]]
        valuesSE = DF_p[DF_p.columns[(DF_p.columns.str.endswith(' down')&~(DF_p.columns.str.startswith('proteomics')))]]
    if mode!='proteomics_only':
        #make sure all CTs in DF and make sure order is correct
        CT_order = get_CT_order(n_clusters,path_filtered_data,True) # endothelial and murine are missing!!
        #CT_order = [str1.replace('-','_').replace(' ','_') for str1 in CT_order]
        #if (len(valuesNW.columns)<len(CT_order)*2) | (len(valuesSE.columns)<len(CT_order)*2):
        #    for ct in CT_order:
        #        if ct+'_up' not in valuesNW.columns.tolist():
        #            valuesNW[ct+'_up'] = np.nan #want them to be grey
        #        if ct+'_down' not in valuesSE.columns.tolist():
        #            valuesSE[ct+'_down'] = np.nan #want them to be grey
        if (len(valuesNW.columns)<len(CT_order)) | (len(valuesSE.columns)<len(CT_order)):
            for ct in CT_order:
                if ct+' up' not in valuesNW.columns.tolist():
                    valuesNW[ct+' up'] = np.nan #want them to be grey
                if ct+' down' not in valuesSE.columns.tolist():
                    valuesSE[ct+' down'] = np.nan #want them to be grey   
    #resort columns such that non-neuronal are always last (before proteomics)
    if any(valuesNW.columns.str.startswith('Astrocytes')) & any(valuesSE.columns.str.startswith('Astrocytes')):
        valuesNW = resort_columns_in_df(valuesNW)
        valuesSE = resort_columns_in_df(valuesSE)
    
    #mask values to plot nan's
    valuesNW_masked = np.ma.masked_invalid(valuesNW)
    valuesSE_masked = np.ma.masked_invalid(valuesSE)

    #mask values to plot significant p-value dots
    valuesNW_masked_sign = valuesNW[valuesNW<=0.05]
    valuesSE_masked_sign = valuesSE[valuesSE<=0.05]

    return [valuesNW, valuesSE], [valuesNW_masked, valuesSE_masked], [valuesNW_masked_sign, valuesSE_masked_sign]

def get_unique_values_in_list_while_preserving_order(input_list):
    _, idx = np.unique(input_list, return_index=True)
    idx_sorted = sorted(idx)
    unique_list_with_preserved_order = [input_list[i] for i in idx_sorted]
    
    return unique_list_with_preserved_order

def resort_columns_in_df(valuesNW):
    #To DO: continue here
    columns = valuesNW.columns
    exc_col = columns[columns.str.startswith('Exc')].tolist()
    inh_col = columns[columns.str.startswith('Inh')].tolist()
    nonneu_col = columns[~columns.str.startswith('Exc') & ~columns.str.startswith('Inh') & ~columns.str.startswith('proteomics')].tolist()
    proteomics_col = columns[columns.str.startswith('proteomics')].tolist()
    new_columns = exc_col + inh_col + nonneu_col + proteomics_col
    if not len(new_columns)==len(columns):
        print("length of new column order differs from length of old columns order. Fix this first!")
    valuesNW = valuesNW.reindex(columns = new_columns)
    return valuesNW

def get_gene_lists_with_qvalues(path_project,version):
    os.chdir(path_project+'data/Gene_lists/')
    if version =='v1' or version=='v2':
        #from genematrix:
        df_gm = pd.read_csv("geneMatrix_"+version+".tsv",sep='\t')
        df_gm.head()
    else:
        genelists=[]
        genelist_names=[]
        q_values=[]
    if version =='v1':
        list_autism_wes, qval_autism_wes = get_list_and_qval(df_gm,'asd.wes.qvalue','gene_name')
        list_scz_wes, qval_scz_wes = get_list_and_qval(df_gm,'scz.wes.qvalue','gene_name')
        list_adhd2018, qval_adhd2018 = get_list_and_qval(df_gm,'adhd2018.gene.P','gene_name')
        list_an2019, qval_an2019 = get_list_and_qval(df_gm,'an2019.gene.P','gene_name')
        list_asd2019, qval_asd2019 = get_list_and_qval(df_gm,'asd2019.gene.P','gene_name')
        list_bip2019, qval_bip2019 = get_list_and_qval(df_gm,'bip2019.gene.P','gene_name')
        list_mdd2019, qval_mdd2019 = get_list_and_qval(df_gm,'mdd2019.gene.P','gene_name')
        list_scz2021, qval_scz2021 = get_list_and_qval(df_gm,'scz2021.gene.P','gene_name')
        list_iq2018, qval_iq2018 = get_list_and_qval(df_gm,'iq2018.gene.P','gene_name')
    
        #from Singh et al --> might have to update this
        df_Singh_SCZ = pd.read_excel('GWAS_scz_Singh2022.xlsx',sheet_name='Table S5 - Gene Results',engine='openpyxl') 
        list_Singh_SCZ, qval_Singh_SCZ = get_list_and_qval(df_Singh_SCZ,'P ca/co (comb), no gnomAD','Gene Symbol')
        
        #from Gandal et al
        df_Gandal_BD = pd.read_csv('Gandal_sign_DEGS_BD.csv')
        list_Gandal_BD, qval_Gandal_BD = get_list_and_qval(df_Gandal_BD,'DGE.fdr','gene_name')
        df_Gandal_ASD = pd.read_csv('Gandal_sign_DEGS_ASD.csv')
        list_Gandal_ASD, qval_Gandal_ASD = get_list_and_qval(df_Gandal_ASD,'DGE.fdr','gene_name')
        df_Gandal_SCZ = pd.read_csv('Gandal_sign_DEGS_SCZ.csv')
        list_Gandal_SCZ, qval_Gandal_SCZ = get_list_and_qval(df_Gandal_SCZ,'DGE.fdr','gene_name')
    
        #from Fu et al
        df_Fu_wes = pd.read_csv('wes_Fu2021_autism.tsv',sep='\t')
        list_Fu_autism_wes, qval_Fu_autism_wes = get_list_and_qval(df_Fu_wes,'FDR_TADA_ASD','gene')
        list_Fu_dd_wes, qval_Fu_dd_wes = get_list_and_qval(df_Fu_wes,'FDR_TADA_DD','gene')
        list_Fu_ndd_wes, qval_Fu_ndd_wes = get_list_and_qval(df_Fu_wes,'FDR_TADA_NDD','gene')
        
        genelists = [list_autism_wes,list_scz_wes,list_adhd2018,list_an2019,list_asd2019,list_bip2019,list_mdd2019,list_scz2021,list_iq2018,list_Singh_SCZ,list_Gandal_BD,list_Gandal_ASD,list_Gandal_SCZ,list_Fu_autism_wes,list_Fu_dd_wes,list_Fu_ndd_wes]
        genelist_names = ['list_asd.wes','list_scz.wes','list_adhd2018','list_an2019','list_asd2019','list_bip2019','list_mdd2019','list_scz2021','list_iq2018', 'list_Singh_SCZ','list_Gandal_BD','list_Gandal_ASD', 'list_Gandal_SCZ','list_Fu_autism_wes','list_Fu_dd_wes','list_Fu_ndd_wes']
        q_values = [qval_autism_wes,qval_scz_wes, qval_adhd2018, qval_an2019, qval_asd2019, qval_bip2019, qval_mdd2019, qval_scz2021, qval_iq2018, qval_Singh_SCZ,qval_Gandal_BD,qval_Gandal_ASD,qval_Gandal_SCZ,qval_Fu_autism_wes,qval_Fu_dd_wes,qval_Fu_ndd_wes]
    elif version=='v2':
        list_asd2020_wes, qval_asd2020_wes = get_list_and_qval(df_gm,'asd2020.wes.qvalue','gene_name')#	as above, but q-values not T/F
        list_asd2021_wes, qval_asd2021_wes = get_list_and_qval(df_gm,'asd2021.wes.qvalue','gene_name')#	as above, but q-values not T/F
        list_dd2021_wes, qval_dd2021_wes = get_list_and_qval(df_gm,'dd2021.wes.qvalue','gene_name')# as above, but q-values not T/F
        list_ndd2021_wes, qval_ndd2021_wes = get_list_and_qval(df_gm,'ndd2021.wes.qvalue','gene_name')#	as above, but q-values not T/F
        list_scz2022_wes, qval_scz2022_wes = get_list_and_qval(df_gm,'scz2022.wes.qvalue','gene_name')#	as above, but q-values not T/

        genelists = [list_asd2020_wes,list_asd2021_wes,list_dd2021_wes, list_ndd2021_wes,list_scz2022_wes]
        genelist_names = ['list_asd2020.wes','list_asd2021.wes','list_dd2021.wes','list_ndd2021.wes','list_scz2022.wes']
        q_values = [qval_asd2020_wes,qval_asd2021_wes, qval_dd2021_wes, qval_ndd2021_wes, qval_scz2022_wes]
    
    return genelists, genelist_names, q_values

def get_list_and_qval(df,qval_col_name,gene_col_name):
    q_val_cutoff=0.1
    gene_list = df[df[qval_col_name]<=q_val_cutoff][gene_col_name]
    qval_list = df[df[qval_col_name]<=q_val_cutoff][qval_col_name]
    return gene_list,qval_list

def get_genelist_DEG_DF(DF_DEG, genelist, genelist_name, q_values, string1):
    DF_DEG_c = DF_DEG.copy()
    DF_DEG_c = get_short_CT_names(DF_DEG_c,"celltype")
    #build DF genelist
    if genelist_name=='list_synaptic_genes':
        #merge based on ensgids
        DF_DEG_c[['Gene_short','ensgid']] = DF_DEG_c.Gene.str.split('_',expand=True)
        DF_genelist = pd.DataFrame(data = {'ensgid': genelist})
        DF_genelist["ensgid"] = DF_genelist["ensgid"].str.strip()
        DF_genelist[genelist_name] = 0.0001#ignored later
        DF = pd.merge(DF_DEG_c, DF_genelist, on='ensgid')
    else:
        #merge based on gene name short
        if len(q_values)==0:
            DF_genelist = pd.DataFrame(data = {'Gene_short': genelist})
            DF_genelist[genelist_name] = 0.0001#ignored later
        else:
            DF_genelist = pd.DataFrame(data = {'Gene_short': genelist, genelist_name: q_values})
        #drop nans
        DF_DEG_c.reset_index(inplace=True)
        DF_genelist.dropna(subset = ["Gene_short"], inplace=True)
        DF_genelist.dropna(subset = [genelist_name], inplace=True)
        DF_genelist["Gene_short"] = DF_genelist["Gene_short"].str.strip()
        #DF_DEG = DF_DEG.reset_index(level=["Gene_short"])
        DF = pd.merge(DF_DEG_c, DF_genelist, on='Gene_short')
    #[item for item in np.unique(DF['Gene_short']) if item not in DF_genelist["Gene_short"].tolist()]
    DF.dropna(subset = ['padj'], inplace=True)
    DF.set_index(keys='Gene_short',inplace=True)
    
    DF.rename(columns={genelist_name: string1+genelist_name},inplace=True)
    #neg log10 transformation
    DF['-log10(padj)'] = -np.log10(DF['padj'])
    DF['-log10('+string1+genelist_name+')'] = -np.log10(DF[string1+genelist_name])
    return(DF)

def get_results_DF_for_sign_genes(DF_p_vals_fdr_corr,celltype_folders,alpha_val):
    bool_first = True
    DF_DESeq2_fdr_corr = pd.DataFrame([],columns = ['Gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'padj_fdr', 'Cell type'])
    for k,ctf in enumerate(celltype_folders):
        os.chdir(ctf)
        #get adjusted p-values for all genes tested:
        if os.path.isfile('df_results_no_shrinkage.csv'):
            df_DESeq2_all_genes = pd.read_csv('df_results_no_shrinkage.csv',sep=',',names=["Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"])
            df_DESeq2_all_genes.drop(index=0, inplace=True)
            df_DESeq2_all_genes = df_DESeq2_all_genes.astype({"baseMean":'float',"log2FoldChange":'float',"lfcSE":'float',"stat":'float',"pvalue":'float',"padj":'float'})

        CT_i = "_".join(ctf.split('/')[len(ctf.split('/'))-1].split('_')[1:-1])
        #rename column of genes:
        df_DESeq2_all_genes = df_DESeq2_all_genes.rename(columns={"Unnamed: 0": "Gene"})
        #filter for significant genes (after fdr correction):
        if CT_i in DF_p_vals_fdr_corr.columns.tolist():
            genes_sign_fdr = DF_p_vals_fdr_corr[DF_p_vals_fdr_corr[CT_i]<alpha_val]['Gene'].tolist()
            padj_fdr = DF_p_vals_fdr_corr[DF_p_vals_fdr_corr[CT_i]<alpha_val][CT_i]
            if len(genes_sign_fdr)!=0:
                df_DESeq2_sel = df_DESeq2_all_genes[df_DESeq2_all_genes['Gene'].isin(genes_sign_fdr)].copy()
                df_DESeq2_sel['padj_fdr'] = padj_fdr.tolist()
                df_DESeq2_sel['Cell type'] = CT_i
                if bool_first:
                    DF_DESeq2_fdr_corr = df_DESeq2_sel.copy()
                    bool_first = False
                else:
                    DF_DESeq2_fdr_corr = pd.concat([DF_DESeq2_fdr_corr, df_DESeq2_sel])
    return DF_DESeq2_fdr_corr

def get_short_CT_names(DF_DEGs,ct_column):
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column].str.replace("Inhibitory","Inh")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("Excitatory","Exc")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("Layer_","L")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("Layer","L")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("neurons_","")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace(" neurons","")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("_neurons","")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('Oligodendrocyte_progenitor_cells',"OPCs")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('Oligodendrocyte progenitor cells',"OPCs")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('Microglial_cells',"Microglia")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('Microglial cells',"Microglia")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('mural_cells',"mural")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace('mural cells',"mural")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("_and_","/")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace(" and ","/")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("  "," ")
    DF_DEGs[ct_column+"_short"] = DF_DEGs[ct_column+"_short"].str.replace("L ","L")
    return DF_DEGs

def get_DF_sign(file_names_sign, folder):
    bool_first_sign = True
    DF_sign = []
    for fs in file_names_sign:
        DF = pd.read_csv(folder+fs,sep=';') 
        DF['Celltype'] = fs[:-5]
        if bool_first_sign:
            DF_sign = DF.copy()
            bool_first_sign = False
        else:
            DF_sign = pd.concat([DF_sign, DF])
            
    return(DF_sign)

def get_number_sign_pw_per_gene(DF_sign,directory,options):
    #TO Do: what happens if cts empty?
    #for each gene count how often it occurs in 'genes.in.geneset.tested'
    if options["opt_up_and_down_sep"]==True:
        reg_strings = ['up','down'] 
    else:
        reg_strings = ['up_or_down'] 
    cts = [ s.split('_'+reg_strings[0])[0] for s in DF_sign['TestVar'].unique().tolist() if s.endswith(reg_strings[0])]
    n_genes_reg = [0]*len(reg_strings)
    if len(cts) == 0:
        DF_sign_pw_per_gene=[]
    else:
        for ct in cts:
            N_genes=[]
            genes = []
            regulation = []
            for reg_id,reg_str in enumerate(reg_strings):
                #create list of ensgids that were part of a pathway
                list_ensgid_pws = []
                for list_ct_i in DF_sign[DF_sign['TestVar']==ct+'_'+reg_str]['genes.in.geneset.tested']:
                    list_ensgid_pws_tmp = list_ct_i.split('_')
                    list_ensgid_pws_tmp = [x for x in list_ensgid_pws_tmp if x!='NA']
                    list_ensgid_pws += list_ensgid_pws_tmp
                #create lists of gene and regulation strings
                #get genes:
                file_gene_list = directory + 'gene_list_' + ct + '_' + reg_str + '.csv'
                genes_tmp = pd.read_csv(file_gene_list,header=None,names=['Gene'])['Gene'].tolist()
                genes = genes + genes_tmp
                regulation = regulation + [reg_str]*len(genes_tmp)
                n_genes_reg[reg_id] = len(genes_tmp)

                #count how often a gene was part of a pathway for this ct
                N_genes_tmp = [0]*n_genes_reg[reg_id]
                import operator as op
                for i,deg_item in enumerate(np.array(genes)[np.array(regulation)==reg_str].tolist()):
                    N_genes_tmp[i] = op.countOf(list_ensgid_pws, deg_item)
                N_genes = N_genes+N_genes_tmp

            if ct==cts[0]:
                DF_sign_pw_per_gene = pd.DataFrame(data={"Gene":genes, "regulation":regulation, "Number_significant_pathways_per_gene": N_genes, "celltype":(len(genes))*[ct]}) 
            else:
                DF_tmp = pd.DataFrame(data={"Gene":genes, "regulation":regulation, "Number_significant_pathways_per_gene": N_genes, "celltype":(len(genes))*[ct]}) 
                DF_sign_pw_per_gene = pd.concat([DF_sign_pw_per_gene,DF_tmp])
    return DF_sign_pw_per_gene

def get_full_DEG_results_folder(path_project,n_cl,dm):
    _, path_results = get_paths(path_project,'DEG_visualization')
    if n_cl == 3:
        DEG_results_folder = path_results + str(n_cl)+'_classes/'+'design_with_'+dm +"/"
    else:
        DEG_results_folder = path_results + str(n_cl)+'_CTs/'+'design_with_'+dm +"/"

    return DEG_results_folder

def get_sign_Genes_with_sign_pathways_per_CT(DF_sign,path_project, n_cl, dm,alpha_val_DEGs,opt_up_and_down_sep):
    if opt_up_and_down_sep:
        reg_strings = ['_down','_up']
    else:
        reg_strings = ['']
    #for each DEG in 'genes.in.geneset.tested' extract pathway for significant (in GSA analysis) pathways
    #get DEGs
    DEG_results_folder = get_full_DEG_results_folder(path_project,n_cl,dm)
    
    DF_pws = DF_sign[DF_sign['subgroup'].str.startswith('GO:')]
    #get DEGs
    DF_p_vals = pd.read_csv(DEG_results_folder+'df_p_vals_adj_sum_per_celltype.csv')
    #change '-' in colnames to '_'
    DF_p_vals.columns = DF_p_vals.columns.astype(str).str.replace("-", "_")
    cols = DF_p_vals.columns
    CTs = cols[np.logical_and(np.logical_and(np.logical_and(~cols.str.startswith('bool_'),~cols.str.startswith('Unname')),~cols.str.startswith('Gene')),~cols.str.startswith('ensgid'))]
    bool_first = True
    for ct in CTs:
        DEGs_up=DF_p_vals[np.logical_and(DF_p_vals['bool_'+ct+'_upregulated']==True, DF_p_vals[ct]<=alpha_val_DEGs)]['ensgid'].tolist()
        DEGs_down=DF_p_vals[np.logical_and(DF_p_vals['bool_'+ct+'_upregulated']==False, DF_p_vals[ct]<=alpha_val_DEGs)]['ensgid'].tolist()
        list_ensgid_pws = []
        for reg_str in reg_strings:
            # loop through the rows using iterrows()
            for index, row in DF_pws[DF_pws['TestVar']==ct+reg_str+'_significant'].iterrows():
                list_ensgid_pws_tmp = row['genes.in.geneset.tested'].split('_')
                list_ensgid_pws_tmp = [x for x in list_ensgid_pws_tmp if x!='NA']
                #filter genes for DEGs
                if reg_str == "_up":
                    ensgids = list(set(list_ensgid_pws_tmp).intersection(set(DEGs_up)))
                    reg_str_list = [reg_str[1:]]*len(ensgids)
                elif reg_str == "_down":
                    ensgids = list(set(list_ensgid_pws_tmp).intersection(set(DEGs_down)))
                    reg_str_list = [reg_str[1:]]*len(ensgids)
                else:
                    ensgids = list(set(list_ensgid_pws_tmp).intersection(set(DEGs_down+DEGs_up)))
                    reg_str_list = [reg_str]*len(ensgids)
                pathway = [row['geneset']]*len(ensgids)
                GO_category = [row['subgroup']]*len(ensgids)
                if bool_first:
                    DF = pd.DataFrame(data={"significant Genes":ensgids, "regulation": reg_str_list, "significant_pathways": pathway, "GO category": GO_category, "Celltype":[ct]*len(ensgids)})
                    bool_first = False
                else:
                    DF_tmp = pd.DataFrame(data={"significant Genes":ensgids, "regulation": reg_str_list, "significant_pathways": pathway, "GO category": GO_category, "Celltype":[ct]*len(ensgids)})
                    DF = pd.concat([DF,DF_tmp])
    #To DO:
    #only keep genes in DF that appear in at least two different pathways
    #make significant pathways into two columns: GO ID significant pathway, significant pathway name
    DF[["GO_ID_significant_pathway", "significant_pathway_name"]] = DF.significant_pathways.str.split(" ",expand=True, n=1)
    DF.set_index("significant Genes",inplace=True)
    bool_first = True
    for ct in DF["Celltype"].unique():
        for reg in reg_strings:
            DF_sel = DF[np.logical_and(DF["Celltype"]==ct, DF["regulation"]==reg)]
            DF_filtered_tmp = DF_sel.loc[DF_sel.index.value_counts()!=1]
            if bool_first:
                DF_filtered = DF_filtered_tmp.copy()
                bool_first=False
            else:
                DF_filtered = pd.concat([DF_filtered,DF_filtered_tmp])
    print("entries removed: " + str(len(DF)-len(DF_filtered)) + "("+ str(np.round(100*((len(DF)-len(DF_filtered))/len(DF)))) +"%)")
    return DF_filtered

def get_number_sign_pathways_per_DEG_per_CT(DF_sign,path_project, n_cl, dm,alpha_val_DEGs,opt_filter_for_GOs):

    if opt_filter_for_GOs:
        DF_sign = DF_sign.loc[DF_sign["subgroup"].str.startswith("GO:")]
    #for each DEG count how often it occurs in 'genes.in.geneset.tested'
    
    #get DEGs
    DEG_results_folder = get_full_DEG_results_folder(path_project,n_cl,dm)

    #get DEGs
    DF_p_vals = pd.read_csv(DEG_results_folder+'df_p_vals_adj_sum_per_celltype.csv')
    #change '-' in colnames to '_'
    DF_p_vals.columns = DF_p_vals.columns.astype(str).str.replace("-", "_")
    cols = DF_p_vals.columns
    CTs = cols[np.logical_and(np.logical_and(np.logical_and(~cols.str.startswith('bool_'),~cols.str.startswith('Unname')),~cols.str.startswith('Gene')),~cols.str.startswith('ensgid'))]
    for ct in CTs:
        DEGs_up=DF_p_vals[np.logical_and(DF_p_vals['bool_'+ct+'_upregulated']==True, DF_p_vals[ct]<=alpha_val_DEGs)]['ensgid'].tolist()
        DEGs_down=DF_p_vals[np.logical_and(DF_p_vals['bool_'+ct+'_upregulated']==False, DF_p_vals[ct]<=alpha_val_DEGs)]['ensgid'].tolist()
        list_ensgid_pws = []
        for reg_str in ['down','up']:
            for list_ct_i in DF_sign[DF_sign['TestVar']==ct+'_'+reg_str+'_significant']['genes.in.geneset.tested']:
                list_ensgid_pws_tmp = list_ct_i.split('_')
                list_ensgid_pws_tmp = [x for x in list_ensgid_pws_tmp if x!='NA']
                list_ensgid_pws += list_ensgid_pws_tmp
        N_DEGs_up = [0]*len(DEGs_up)
        N_DEGs_down = [0]*len(DEGs_down)
        import operator as op
        for i,deg_item in enumerate(DEGs_up):
            N_DEGs_up[i] = op.countOf(list_ensgid_pws, deg_item)
        for i,deg_item in enumerate(DEGs_down):
            N_DEGs_down[i] = op.countOf(list_ensgid_pws, deg_item)
        if ct==CTs[0]:
            DF_sign_pw_per_DEG = pd.DataFrame(data={"DEGs":DEGs_up+DEGs_down, 
                                                    "regulation":len(DEGs_up)*['up']+len(DEGs_down)*['down'], 
                                                    "Number_significant_pathways_per_gene": N_DEGs_up+N_DEGs_down, 
                                                    "celltype":(len(DEGs_up)+len(DEGs_down))*[ct]
                                                    }) 
        else:
            DF_tmp = pd.DataFrame(data={"DEGs":DEGs_up+DEGs_down, 
                                        "regulation":len(DEGs_up)*['up']+len(DEGs_down)*['down'], 
                                        "Number_significant_pathways_per_gene": N_DEGs_up+N_DEGs_down, 
                                        "celltype":(len(DEGs_up)+len(DEGs_down))*[ct]
                                        }) 
            DF_sign_pw_per_DEG = pd.concat([DF_sign_pw_per_DEG,DF_tmp])

    return DF_sign_pw_per_DEG

def get_statistics_pathways(DF_up, DF_down, df_number_pws_not_tested, p_val_cutoff,opt_proteomics_up_and_down_sep):
    #TO DO: use df_number_pws_tested
    df_number_pws_not_tested_t = df_number_pws_not_tested.set_index('celltype_name').transpose()
    #extract statistics for number of sign pathways etc
    #nans = pd.concat([DF_up.isna().sum(),DF_down.isna().sum()])
    groups = np.unique(DF_up.index.get_level_values('group').tolist())
    if len(groups)>1: #take total
        n_pw_not_tested = df_number_pws_not_tested_t.loc['total']
    else:
        n_pw_not_tested = df_number_pws_not_tested_t.loc[groups[0]]
        
    #n_pw_tested = np.shape(DF_up)[0]
    n_sign = pd.concat([(DF_up<=p_val_cutoff).sum(),(DF_down<=p_val_cutoff).sum()])
    n_nsign = pd.concat([(DF_up>p_val_cutoff).sum(),(DF_down>p_val_cutoff).sum()])
    DF_n_pw = pd.concat({'significant':n_sign,'non-significant':n_nsign,'not tested': n_pw_not_tested},axis=1, ignore_index=False)
    
    if opt_proteomics_up_and_down_sep==False:
        DF_n_pw.drop(index=['proteomics_up'],inplace=True)
        DF_n_pw.rename(index={'proteomics_down':'proteomics'},inplace=True)
    DF_n_pw = DF_n_pw.rename_axis('celltype_name').reset_index(level=0)

    DF_n_pw[['Celltype', 'regulation']] = DF_n_pw["celltype_name"].str.rsplit(pat=r'_',n=1,expand=True)
    DF_n_pw.loc[DF_n_pw['Celltype']=='proteomics','regulation']='deregulated'
    DF_n_pw.drop(columns='celltype_name',inplace=True)
    #make CT index:
    #DF_n_pw = DF_n_pw.set_index('Celltype')
    DF_n_pw['tested'] = DF_n_pw['significant']+DF_n_pw['non-significant']
    DF_n_pw['Percentage significant'] = 100*DF_n_pw['significant']/DF_n_pw['tested']
    DF_n_pw['Percentage tested'] = 100*DF_n_pw['tested']/(DF_n_pw['tested'] + DF_n_pw['not tested'])
    DF_n_pw['Percentage significant'][DF_n_pw['tested'] == 0] = 0
    
    #remove some dashes
    DF_n_pw["Celltype"] = [ct.replace('2_3','2-3').replace('5_6','5-6').replace('3_6','3-6').replace('3_4','3-4').replace('_',' ') for ct in DF_n_pw["Celltype"].to_list()]
    #get short cell type names
    DF_n_pw = get_short_CT_names(DF_n_pw,"Celltype")

    return DF_n_pw


def get_statistics_pathways_all(DF_all, df_number_pws_not_tested, p_val_cutoff, pval_col_to_extract):
    # to do: % sign out of pws tested, all GOs together
    df_number_pws_not_tested_t = df_number_pws_not_tested.set_index('celltype_name').transpose()
    #extract statistics for number of sign pathways etc
    #nans = pd.concat([DF_up.isna().sum(),DF_down.isna().sum()])
    groups = np.unique(DF_all.index.get_level_values('group').tolist())
    if len(groups)>1: #take total
        n_pw_not_tested = df_number_pws_not_tested_t.loc['total']
    else:
        n_pw_not_tested = df_number_pws_not_tested_t.loc[groups[0]]
        
    #n_pw_tested = np.shape(DF_up)[0]
    n_sign = (DF_all<=p_val_cutoff).sum()
    n_nsign = (DF_all>p_val_cutoff).sum()
    DF_n_pw = pd.concat({'significant':n_sign,'non-significant':n_nsign,'not tested': n_pw_not_tested},axis=1, ignore_index=False)
    DF_n_pw = DF_n_pw.rename_axis('celltype_name').reset_index(level=0)

    DF_n_pw['regulation']='deregulated'
    DF_n_pw.rename(columns={'celltype_name':'Celltype'},inplace=True)
    #DF_n_pw.drop(columns='celltype_name',inplace=True)
    #make CT index:
    #DF_n_pw = DF_n_pw.set_index('Celltype')
    DF_n_pw['tested'] = DF_n_pw['significant']+DF_n_pw['non-significant']
    DF_n_pw['Percentage significant'] = 100*DF_n_pw['significant']/DF_n_pw['tested']
    DF_n_pw['Percentage tested'] = 100*DF_n_pw['tested']/(DF_n_pw['tested'] + DF_n_pw['not tested'])
    DF_n_pw['Percentage significant'][DF_n_pw['tested'] == 0] = 0
    
    #reove some dashes
    DF_n_pw["Celltype"] = [ct.replace('2_3','2-3').replace('5_6','5-6').replace('3_6','3-6').replace('3_4','3-4').replace('_',' ') for ct in DF_n_pw["Celltype"].to_list()]
    #get short cell type names
    DF_n_pw = get_short_CT_names(DF_n_pw,"Celltype")

    return DF_n_pw

def get_PCs_for_pathways(DF_sign,variable_gene_names,variables_groups):
    #PCA of pathways:
    #reshape data:
    DF_sign.reset_index(inplace=True,drop=True)
    n_pathways = np.shape(DF_sign)[0]
    genes=[]
    for p in range(0,n_pathways):
        genes = genes + DF_sign.iloc[p][variable_gene_names].split('_')
    genes_list = np.unique(genes).tolist()
    genes_list = [x for x in genes_list if x != 'NA']
    df = pd.DataFrame(columns=genes_list,index=DF_sign.index.tolist())
    for r in df.index.tolist():
        genes_r = DF_sign.iloc[r][variable_gene_names].split('_')
        genes_r = [x for x in genes_r if x != 'NA']
        for g1 in genes_r:
            df.at[r,g1]=1
        for g0 in list(set(genes_list)-set(genes_r)):
            df.at[r,g0]=0
        #df.at[r,genes_r]=1
        #df.at[r,list(set(genes_list)-set(genes_r))]=0
    #add go, pathway, ct info
    df[variables_groups[0]] = DF_sign[variables_groups[0]]
    df[variables_groups[1]] = DF_sign[variables_groups[1]]
    df['Celltype'] = DF_sign['Celltype']
    df.set_index([variables_groups[0],variables_groups[1],'Celltype'],inplace=True)
    #visualize:
    
    #perform PCA
    features = df.columns.tolist()
    # Separating out the features
    x = df.loc[:, features].values
    # Standardizing the features
    x = StandardScaler().fit_transform(x)
    
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'],index=df.index)
    
    return(principalDf)

def get_DFs_all(file_names_all, pval_col_to_extract, alpha_val_DEGs, folder_transcriptomics, path_project, n_cl, dm, opt_proteomics, results_path_proteomics,path_results,opt_load_dataframe,opt_filter_stats_for_GOs):
    #TO DO: add percentage up-regulated DEGs in each pathway to later determine plainly up-reg pathways, plainly down-reg pathways, mixed pathways
    if opt_load_dataframe==False:
        #proteomics is never seperately tested in this case
        bool_first_all = True
        
        DEG_results_folder = get_full_DEG_results_folder(path_project,n_cl,dm)
    
        #get DEGs
        DF_p_vals = pd.read_csv(DEG_results_folder+'df_p_vals_adj_sum_per_celltype.csv')
        #change '-' in colnames to '_'
        DF_p_vals.columns = DF_p_vals.columns.astype(str).str.replace("-", "_")
        
        for fa in file_names_all:
            if fa.endswith('__all.csv'):
                #for p-val heatmap:
                #make sure nothing gets stored twice:
                df_tmp = pd.DataFrame()
                df_tra = pd.DataFrame()
                
                df_tmp = pd.read_csv(folder_transcriptomics+fa,sep=';')
                df_tra = df_tmp[['TestVar','group','subgroup','geneset',pval_col_to_extract,'genes.in.geneset.tested']].copy()
                
                if "Patient_clustering/" in folder_transcriptomics:
                    #get cell type names:
                    df_tra['celltype_name'] = df_tra['TestVar']
                    #get list of genes for current cell type:
                    ct = df_tra['celltype_name'].unique().tolist()[0]
                    #if opt_up_and_down_separately==True:
                    #needs to be a list, not dataframe 
                    genes_patient_clustering = pd.read_csv(folder_transcriptomics+'gene_list_'+ct+'.csv',header=None)[0].tolist()
                    #determine number of up-regulated and number of downregulated genes in pathways
                    if '__down_all' in fa:
                        reg_string = 'down'
                    elif '__up_all' in fa:
                        reg_string = 'up'
                    df_tra['number '+reg_string+'regulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str = reg_string, alpha_val_DEGs=[], DF_p_vals=[], gene_list = genes_patient_clustering, opt_up_and_down_separately=[], axis=1)
                    #else:
                        #genes_patient_clustering = pd.read_csv(folder_transcriptomics+'gene_list_'+ct+'up_or_down_all.csv')
                else:
                    df_tra['celltype_name'] = df_tra['TestVar'].str.split('_significant',n=1).str[0]
                    #determine number of up-regulated and number of downregulated DEGs in pathways
                    df_tra['number downregulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str = 'down', alpha_val_DEGs=alpha_val_DEGs, DF_p_vals=DF_p_vals, gene_list =[], opt_up_and_down_separately=False, axis=1)
                    df_tra['number upregulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str='up', alpha_val_DEGs=alpha_val_DEGs, DF_p_vals=DF_p_vals, gene_list = [], opt_up_and_down_separately=False, axis=1)
                
                # track number of pathways tested per CT:
                if opt_filter_stats_for_GOs:
                    pw_categories = ["gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF"]
                    n_pws_total = [7658,1006,1738]
                else:
                    pw_categories = ["gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF",
                                    'disease', 'dna-rna-binding', 'drugs', 'evo-constraint', 'expression', 'gene-annotation',
                                    "gene-ontology_HPO","gene-ontology_mitochondria",
                                    'other', 'synaptic-gene-ontology', 'wikipathways']
                    n_pws_total = [7658,1006,1738,10,707,12,6,312,6, 5071,1, 4,293,692]
                    
                # pw_categories = ['disease', 'dna-rna-binding', 'drugs', 'evo-constraint', 'expression', 'gene-annotation',
                #                  "gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF","gene-ontology_HPO","gene-ontology_mitochondria",
                #                  'other', 'synaptic-gene-ontology', 'wikipathways']
                # n_pws_total = [10,707,12,6,312,6, 7658,1006,1738,5071,1, 4,293,692]
                if (bool_first_all):
                    df_number_pws_not_tested = pd.DataFrame(data={'celltype_name':np.unique(df_tra['celltype_name']).tolist()[0], 'total': np.sum(n_pws_total)-len(df_tra)},index=[0])
                    for i,gr in enumerate(pw_categories):
                        df_number_pws_not_tested[gr] = n_pws_total[i]-sum(df_tra['group']==gr)
                else:
                    data=[np.unique(df_tra['celltype_name']).tolist()[0],np.sum(n_pws_total)-len(df_tra)]
                    for i,gr in enumerate(pw_categories):
                        data.append(n_pws_total[i]-sum(df_tra['group']==gr))
                    df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                    
                #print('number of duplicates for '+fa)
                #print(np.shape(df_tra)[0]-np.shape(df_tra.drop_duplicates(keep='first',inplace = False, ignore_index=False))[0])
                if bool_first_all:
                    DF_tra = df_tra.copy()
                    bool_first_all = False
                else:
                    DF_tra = pd.concat([DF_tra, df_tra.copy()])
            else:
                continue
        
        #DF_tra_down.drop_duplicates(keep='first',inplace = True, subset=index_cols+['celltype_name'], ignore_index=True)
        #DF_tra_up.drop_duplicates(keep='first',inplace = True, subset=index_cols+['celltype_name'], ignore_index=True)
        if opt_proteomics==True and os.path.isdir(results_path_proteomics):
            #add proteomics result
            # get all file names:
            bool_proteomics_first = True
            files = os.listdir(results_path_proteomics)
            for f in files:
                if not f.endswith('both_all.csv'):
                    continue
                df_tmp = pd.read_csv(results_path_proteomics+f,sep=';')
                df_pro = df_tmp[['group','subgroup','geneset',pval_col_to_extract]]
                if bool_proteomics_first == True:
                    df_pro['celltype_name']='proteomics'
                    #put at the end of trans up and down:
                    DF_all = pd.concat([DF_tra,df_pro.copy()])
                    bool_proteomics_first = False
                else:
                    DF_all = pd.concat([DF_all,df_pro.copy()])
                #add number of significant pathways:
                data=[np.unique(DF_all['celltype_name'][DF_all['celltype_name'].str.startswith('prot')]).tolist()[0],np.sum(n_pws_total)-len(df_pro)]
                for i,gr in enumerate(pw_categories):
                    data.append(n_pws_total[i]-sum(df_pro['group']==gr))
                df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
        else:
            DF_all = DF_tra.copy()
            
        # after new cell type clustering this should no longer be necessary
        #DF_all_up = DF_all_up.drop_duplicates(keep='first',subset=index_cols+['celltype_name'],ignore_index=True)
        #DF_all_down = DF_all_down.drop_duplicates(keep='first',subset=index_cols+['celltype_name'], ignore_index=True)
        if opt_filter_stats_for_GOs:
            DF_all = DF_all.loc[DF_all["subgroup"].str.startswith("GO:")]

        #reindex:
        DF_all.reset_index(drop=True,inplace=True)
        
        DF_all.to_csv(path_results+'DF_all.csv')
        df_number_pws_not_tested.to_csv(path_results+'df_number_pws_not_tested.csv')
    else:
        DF_all=pd.read_csv(path_results+'DF_all.csv')
        df_number_pws_not_tested=pd.read_csv(path_results+'df_number_pws_not_tested.csv')
    return DF_all, df_number_pws_not_tested

def get_number_genes_per_pathway(row,reg_str, alpha_val_DEGs, DF_p_vals, gene_list, opt_up_and_down_separately):
    #called from get_DFs_all()
    list_ensgid = row['genes.in.geneset.tested'].split('_')
    list_ensgid = [x for x in list_ensgid if x!='NA']
    if len(gene_list)==0:
        #select from DF_p_vals all ensgids in list_ensgid
        DF_p_vals_sel = DF_p_vals[DF_p_vals['ensgid'].isin(list_ensgid)]
        #check how many are significant and up/down/in total
        if opt_up_and_down_separately:
        #pick cell type column of interest:
            if 'bool_'+row['celltype_name']+'regulated' in DF_p_vals_sel.columns:
                n=np.sum(np.logical_and(DF_p_vals_sel['bool_'+row['celltype_name']+'regulated']==False, DF_p_vals_sel[row['celltype_name'].split('_up')[0].split('_down')[0]]<=alpha_val_DEGs))
            else:
                n=0
        else:
            if reg_str == 'up':
                n=np.sum(np.logical_and(DF_p_vals_sel['bool_'+row['celltype_name']+'_upregulated']==True, DF_p_vals_sel[row['celltype_name']]<=alpha_val_DEGs))
            elif reg_str =='down':
                n=np.sum(np.logical_and(DF_p_vals_sel['bool_'+row['celltype_name']+'_downregulated']==False, DF_p_vals_sel[row['celltype_name']]<=alpha_val_DEGs))
            else: #other
                n=np.nan
    else:
        #to do: make sure gene_list is list, no df
        genes_mapping_to_pathway = [gene for gene in gene_list if gene in list_ensgid]
        n=len(genes_mapping_to_pathway)
    return n

def clustermap_pathway_order(df_tmp):
    g = sns.clustermap(df_tmp, yticklabels=1,xticklabels=1,col_cluster = False)
    pathway_list_ordered = df_tmp.index[g.dendrogram_row.reordered_ind].tolist()
    
    return pathway_list_ordered

def get_DFs_all_up_down(file_names_all, pval_col_to_extract, alpha_val_DEGs, folder_transcriptomics, path_project, n_cl, dm, opt_proteomics, opt_proteomics_up_and_down_sep, results_path_proteomics,path_results,opt_load_dataframe,opt_filter_stats_for_GOs):
    if opt_filter_stats_for_GOs:
        add_str = "_GOs_only"
    else:
        add_str = ""
    
    if opt_load_dataframe==False:
        #triangular heatmap:
        #get transcriptomic result:
        bool_first_all_up = True
        bool_first_all_down = True
        
        if "Patient_clustering/" not in folder_transcriptomics:
            #get DEGs for up and down and all cell types
            DEG_results_folder = get_full_DEG_results_folder(path_project,n_cl,dm)
            
            DF_p_vals = pd.read_csv(DEG_results_folder+'df_p_vals_adj_sum_per_celltype.csv')
            #change '-' in colnames to '_'
            DF_p_vals.columns = DF_p_vals.columns.astype(str).str.replace("-", "_")
        
        for fa in file_names_all:
            if fa.endswith('_up_all.csv'):
                reg_str_tra_i='up'
            elif fa.endswith('_down_all.csv'):
                if fa.endswith('_up_or_down_all.csv'): # ignore those files 
                    continue
                else:
                    reg_str_tra_i='down'
            else:
                continue
            #for p-val heatmap:
                
            #make sure nothing gets stored twice:
            df_tmp = pd.DataFrame()
            df_tra = pd.DataFrame()
            
            df_tmp = pd.read_csv(folder_transcriptomics+fa,sep=';')
            df_tra = df_tmp[['TestVar','group','subgroup','geneset','genes.in.geneset.tested',pval_col_to_extract]].copy()  

            if "Patient_clustering/" in folder_transcriptomics:
                #DF_down contains too many columns (also up_or_down.csv files are added) --> fix this!
                df_tra['celltype_name'] = df_tra['TestVar']
                #get list of genes for current cell type:
                ct = df_tra['celltype_name'].unique().tolist()[0]
                genes_patient_clustering = pd.read_csv(folder_transcriptomics+'gene_list_'+ct+'.csv',header=None)[0].tolist()
                #determine number of up-regulated and number of downregulated genes in pathways
                if '__down_all' in fa:
                    reg_string = 'down'
                elif '__up_all' in fa:
                    reg_string = 'up'
                df_tra['number '+reg_string+'regulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str = reg_string, alpha_val_DEGs=[], DF_p_vals=[], gene_list = genes_patient_clustering, opt_up_and_down_separately=[], axis=1)
       
            else:
                df_tra['celltype_name'] = df_tra['TestVar'].str.split('_'+reg_str_tra_i+'_significant',n=1).str[0]
                #add info if up or down and restructure:
                df_tra['celltype_name'] = df_tra['celltype_name']+'_'+reg_str_tra_i
                #determine number of up-regulated and number of downregulated DEGs in pathways
                df_tra['number downregulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str = 'down', alpha_val_DEGs=alpha_val_DEGs, DF_p_vals=DF_p_vals, gene_list = [], opt_up_and_down_separately=True, axis=1)
                df_tra['number upregulated genes']=df_tra.apply(get_number_genes_per_pathway,reg_str = 'up', alpha_val_DEGs=alpha_val_DEGs, DF_p_vals=DF_p_vals, gene_list = [], opt_up_and_down_separately=True, axis=1)

            # track numbe rof pathways tested per CT:
            if opt_filter_stats_for_GOs:
                pw_categories = ["gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF"]
                n_pws_total = [7658,1006,1738]
            else:
                pw_categories = ["gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF",
                                'disease', 'dna-rna-binding', 'drugs', 'evo-constraint', 'expression', 'gene-annotation',
                                "gene-ontology_HPO","gene-ontology_mitochondria",
                                'other', 'synaptic-gene-ontology', 'wikipathways']
                n_pws_total = [7658,1006,1738,10,707,12,6,312,6, 5071,1, 4,293,692]
            # pw_categories = ['disease', 'dna-rna-binding', 'drugs', 'evo-constraint', 'expression', 'gene-annotation',
            #                  "gene-ontology_GO:BP","gene-ontology_GO:CC","gene-ontology_GO:MF","gene-ontology_HPO","gene-ontology_mitochondria",
            #                  'other', 'synaptic-gene-ontology', 'wikipathways']
            # n_pws_total = [10,707,12,6,312,6, 7658,1006,1738,5071,1, 4,293,692]
            if (bool_first_all_down & bool_first_all_up):
                df_number_pws_not_tested = pd.DataFrame(data={'celltype_name':np.unique(df_tra['celltype_name']).tolist()[0], 'total': np.sum(n_pws_total)-len(df_tra)},index=[0])
                for i,gr in enumerate(pw_categories):
                    df_number_pws_not_tested[gr] = n_pws_total[i]-sum(df_tra['group']==gr)
            else:
                data=[np.unique(df_tra['celltype_name']).tolist()[0],np.sum(n_pws_total)-len(df_tra)]
                for i,gr in enumerate(pw_categories):
                    data.append(n_pws_total[i]-sum(df_tra['group']==gr))
                df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                
            #print('number of duplicates for '+fa)
            #print(np.shape(df_tra)[0]-np.shape(df_tra.drop_duplicates(keep='first',inplace = False, ignore_index=False))[0])
            if reg_str_tra_i=='up':
                if bool_first_all_up:
                    DF_tra_up = df_tra.copy()
                    bool_first_all_up = False
                else:
                    DF_tra_up = pd.concat([DF_tra_up, df_tra.copy()])
            elif reg_str_tra_i=='down':
                if bool_first_all_down:
                    DF_tra_down = df_tra.copy()
                    bool_first_all_down = False
                else:
                    DF_tra_down = pd.concat([DF_tra_down, df_tra.copy()])
        
        #DF_tra_down.drop_duplicates(keep='first',inplace = True, subset=index_cols+['celltype_name'], ignore_index=True)
        #DF_tra_up.drop_duplicates(keep='first',inplace = True, subset=index_cols+['celltype_name'], ignore_index=True)
        if opt_proteomics==True and os.path.isdir(results_path_proteomics):
            #add proteomics result
            # get all file names:
            bool_proteomics_up_first = True
            bool_proteomics_down_first = True
            files = os.listdir(results_path_proteomics)
            for f in files:
                if opt_proteomics_up_and_down_sep:
                    if f.endswith('up_all.csv'):
                        reg_str_pro_i = 'up'
                    elif f.endswith('down_all.csv'):
                        reg_str_pro_i = 'down'
                    else:
                        continue
                else:
                    if f.endswith('both_all.csv'):
                        reg_str_pro_i = 'both'
                    else:
                        continue
                        
                df_tmp = pd.read_csv(results_path_proteomics+f,sep=';')
                df_pro = df_tmp[['group','subgroup','geneset',pval_col_to_extract]]
                
                if opt_proteomics_up_and_down_sep:
                    df_pro['celltype_name']='proteomics_'+reg_str_pro_i
                    #directly concat up- and down-reg data frames respectively:
                    if reg_str_pro_i=='up':
                        #put at the end of trans up
                        if bool_proteomics_up_first == True:
                            DF_all_up = pd.concat([DF_tra_up,df_pro.copy()])
                            bool_proteomics_up_first = False
                        else:
                            DF_all_up = pd.concat([DF_all_up,df_pro.copy()])
                        #add number of significant pathways:
                        data=[np.unique(DF_all_up['celltype_name'][DF_all_up['celltype_name'].str.startswith('prot')]).tolist()[0],np.sum(n_pws_total)-len(df_pro)]
                        for i,gr in enumerate(pw_categories):
                            data.append(n_pws_total[i]-sum(df_pro['group']==gr))
                        df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                    elif reg_str_pro_i=='down':
                        #put at the end of trans down
                        if bool_proteomics_down_first == True:
                            DF_all_down = pd.concat([DF_tra_down,df_pro.copy()])
                            bool_proteomics_down_first = False
                        else:
                            DF_all_down = pd.concat([DF_all_down,df_pro.copy()])
                        #add number of significant pathways:
                        data=[np.unique(DF_all_down['celltype_name'][DF_all_down['celltype_name'].str.startswith('prot')]).tolist()[0],np.sum(n_pws_total)-len(df_pro)]
                        for i,gr in enumerate(pw_categories):
                            data.append(n_pws_total[i]-sum(df_pro['group']==gr))
                        df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                else:
                    if bool_proteomics_up_first == True:
                        df_pro['celltype_name']='proteomics_up'
                        #put at the end of trans up and down:
                        DF_all_up = pd.concat([DF_tra_up,df_pro.copy()])
                        bool_proteomics_up_first = False
                    else:
                        DF_all_up = pd.concat([DF_all_up,df_pro.copy()])
                    #add number of significant pathways:
                    data=[np.unique(DF_all_up['celltype_name'][DF_all_up['celltype_name'].str.startswith('prot')]).tolist()[0],np.sum(n_pws_total)-len(df_pro)]
                    for i,gr in enumerate(pw_categories):
                        data.append(n_pws_total[i]-sum(df_pro['group']==gr))
                    df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                    
                    if bool_proteomics_down_first == True:
                        df_pro['celltype_name']='proteomics_down'
                        DF_all_down = pd.concat([DF_tra_down,df_pro.copy()])
                        bool_proteomics_down_first = False
                    else:
                        DF_all_down = pd.concat([DF_all_down,df_pro.copy()])
                    #add number of significant pathways:
                    data=[np.unique(DF_all_down['celltype_name'][DF_all_down['celltype_name'].str.startswith('prot')]).tolist()[0],np.sum(n_pws_total)-len(df_pro)]
                    for i,gr in enumerate(pw_categories):
                        data.append(n_pws_total[i]-sum(df_pro['group']==gr))
                    df_number_pws_not_tested.loc[len(df_number_pws_not_tested.index)] = data
                
            #make sure dfs contain no duplicate rows:

        else:
            DF_all_up = DF_tra_up.copy()
            DF_all_down = DF_tra_down.copy()
            
        # after new cell type clustering this should no longer be necessary
        #DF_all_up = DF_all_up.drop_duplicates(keep='first',subset=index_cols+['celltype_name'],ignore_index=True)
        #DF_all_down = DF_all_down.drop_duplicates(keep='first',subset=index_cols+['celltype_name'], ignore_index=True)
        
        if opt_filter_stats_for_GOs:
            DF_all_up = DF_all_up.loc[DF_all_up["subgroup"].str.startswith("GO:")]
            DF_all_down = DF_all_down.loc[DF_all_down["subgroup"].str.startswith("GO:")]

        #reindex:
        DF_all_up.reset_index(drop=True,inplace=True)
        DF_all_down.reset_index(drop=True,inplace=True)
                
        #save Data frames:
        DF_all_up.to_csv(path_results+'DF_all_up'+add_str+'.csv')
        DF_all_down.to_csv(path_results+'DF_all_down'+add_str+'.csv')
        df_number_pws_not_tested.to_csv(path_results+'df_number_pws_not_tested'+add_str+'.csv')

    else:
        DF_all_up=pd.read_csv(path_results+'DF_all_up'+add_str+'.csv')
        DF_all_down=pd.read_csv(path_results+'DF_all_down'+add_str+'.csv')
        df_number_pws_not_tested=pd.read_csv(path_results+'df_number_pws_not_tested'+add_str+'.csv')
        
    return DF_all_up, DF_all_down, df_number_pws_not_tested

def get_DFs_as_pivot(DF_up, DF_down,pval_col_to_extract,index_cols):
    #make cell type name new columns and pvalues values:
    DF_up = DF_up.pivot(index=index_cols,columns='celltype_name',values=pval_col_to_extract)
    DF_down = DF_down.pivot(index=index_cols,columns='celltype_name',values=pval_col_to_extract)
    #replace komma with dot and make double:
    DF_up=DF_up.astype('str').replace(',','.',regex=True).astype('float')
    DF_down=DF_down.astype('str').replace(',','.',regex=True).astype('float')
    return DF_up, DF_down

def get_DFs_as_pivot_all(DF_all, pval_col_to_extract, index_cols):
    #make cell type name new columns and pvalues values:
    DF_all = DF_all.pivot(index=index_cols,columns='celltype_name',values=pval_col_to_extract)
    #replace komma with dot and make double:
    DF_all=DF_all.astype('str').replace(',','.',regex=True).astype('float')
    return DF_all

def hex_to_rgb(hx, hsl=False):
    #from https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python
    """Converts a HEX code into RGB or HSL.
    Args:
        hx (str): Takes both short as well as long HEX codes.
        hsl (bool): Converts the given HEX code into HSL value if True.
    Return:
        Tuple of length 3 consisting of either int or float values.
    Raise:
        ValueError: If given value is not a valid HEX code."""
    if re.compile(r'#[a-fA-F0-9]{3}(?:[a-fA-F0-9]{3})?$').match(hx):
        div = 255.0 if hsl else 0
        if len(hx) <= 4:
            return tuple(int(hx[i]*2, 16) / div if div else
                         int(hx[i]*2, 16) for i in (1, 2, 3))
        return tuple(int(hx[i:i+2], 16) / div if div else
                     int(hx[i:i+2], 16) for i in (1, 3, 5))
    raise ValueError(f'"{hx}" is not a valid HEX code.')

# Return one 24-bit color value 
def rgb_to_dec(rgb_tuple):
    r = rgb_tuple[0]/255.0
    g = rgb_tuple[1]/255.0
    b = rgb_tuple[2]/255.0
    return (r,b,g)
    #https://stackoverflow.com/questions/25404998/how-can-i-convert-an-rgb-input-into-a-decimal-color-code
    #return (r << 16) + (g << 8) + b

def get_continuous_cmap(hex_list, float_list=None):
    #from: https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def generate_GSA_output_and_statistics(path_info, dm, options, variable_gene_names, n_cl,alpha_val_DEGs,variables_groups,variable_n_genes,background_data_proteomics,directory,j,bg_str,results_path_proteomics):
    #To Do: certain functions of this method should only be called for DEGs
    #for other modes dm empty string and alpha_val_DEGs=[]
    if options["opt_up_and_down_sep"]:
        add_str = 'all_separated'
    else:
        add_str = 'all_together'

    #if path_info["path_results_subdirectory"] exists:
    if os.path.isdir(path_info["path_results_subdirectory"]):
        files = os.listdir(path_info["path_results_subdirectory"])
        if options["opt_up_and_down_sep"]:
            file_names_all = [s for s in files if (s.endswith('_up_all.csv') or s.endswith('_down_all.csv'))]
            file_names_sign = [s for s in files if (s.endswith('_up_sign.csv') or s.endswith('_down_sign.csv'))]
        else:
            file_names_all = [s for s in files if s.endswith('__all.csv')]
            file_names_sign = [s for s in files if s.endswith('__sign.csv')]
        if len(file_names_sign)>0:
            #get data frame of significant pathways listed for various cell types:
            DF_sign = get_DF_sign(file_names_sign, path_info["path_results_subdirectory"])
            groups = np.unique(DF_sign['group'])
            if "Patient_clustering" in directory: #GSA of patient clustering gene tensors
                #TO DO: implement
                #plot significant pathways per CT per gene
                #for each gene count how often it occurs in 'genes.in.geneset.tested'
                
                DF_sign_pw_per_gene = get_number_sign_pw_per_gene(DF_sign,path_info["path_results_subdirectory"],options,False)
                #plot_number_pathways_per_gene(DF_sign_pw_per_gene, options["opt_save"], path_info["path_results_subdirectory"],n_cl,path_info["path_colors"],  "")
                for group_name in groups:
                    DF_sign_group_i = DF_sign[DF_sign['group']==group_name].copy()
                    DF_sign_pw_per_gene_gr_i = get_number_sign_pw_per_gene(DF_sign_group_i,path_info["path_results_subdirectory"],options,False)
                    if len(DF_sign_pw_per_gene_gr_i)>0:
                        plot_number_pathways_per_gene(DF_sign_pw_per_gene_gr_i, options["opt_save"], path_info["path_results_subdirectory"],n_cl,path_info["path_colors"], group_name)
                #UMAP/ PCA taking pseudobulk expression data of gene list genes for each factor and each CT separately. if whole data set: color code for SCZ and CTRL samples
                cts = [s.split('_up')[0] for s in DF_sign['TestVar'].unique().tolist() if s.endswith('up')]
                factor = path_info["path_results_subdirectory"].split('/')[-2]
                plot_PCA_of_pseudobulk_data_for_gene_selection(cts,options, path_info, n_cl, factor, "input_GSA_of_factors","")
                
            else: #GSA of DEGs
                # which biol. functions are using the same disrupted genes?
                # investigate how many DEGs map to two or more pathways with different parent terms: what are the parent terms?
                DF_pw_sharing_DEGs = get_sign_Genes_with_sign_pathways_per_CT(DF_sign,path_info["path_project"], n_cl, dm,alpha_val_DEGs, options["opt_up_and_down_sep"])
                #later with rrvgo info: identify which genes are shared across parent Terms and which parent terms these are: is there a pattern?
                GO_term_str = ['BP','MF','CC']
                for go_str in GO_term_str:
                    #get GO term clustering:
                    rrvgo_info_strings = ['up','down']
                    for rrvgo_info_str in rrvgo_info_strings:
                        rrvgo_filename = 'module_list_rrvgo_all_'+go_str+'_size_'+rrvgo_info_str+'.csv'
                        if os.path.isfile(path_info["path_results_subdirectory"]+rrvgo_filename): 
                            module_info_tmp = pd.read_csv(path_info["path_results_subdirectory"]+rrvgo_filename)
                            module_info_tmp.drop(columns=[c for c in module_info_tmp.columns if c not in ['go','term',"parent","parentTerm"]])
                            if rrvgo_info_str=='up':
                                module_info = module_info_tmp
                            else:
                                module_info = pd.concat([module_info,module_info_tmp])
                    #integrate this info for current GO term in DF_pw_sharing_DEG

                DF_sign_pw_per_DEG = get_number_sign_pathways_per_DEG_per_CT(DF_sign,path_info["path_project"], n_cl, dm, alpha_val_DEGs,options["opt_filter_stats_for_GOs"])
                if options["opt_filter_stats_for_GOs"]:
                    add_str_2 = "_GOs_only"
                else:
                    add_str_2 = "_all_pws_tested"
                #plot number of significant pathways per down- vs. per up-regulated DEG
                plot_number_pathways_per_gene(DF_sign_pw_per_DEG, options["opt_save"], path_info["path_results_subdirectory"], n_cl, path_info["path_colors"],add_str+add_str_2)
                # for each group separately:
                for group_name in groups:
                    DF_sign_group_i = DF_sign[DF_sign['group']==group_name].copy()
                    DF_sign_pw_per_DEG_gr_i = get_number_sign_pathways_per_DEG_per_CT(DF_sign_group_i,path_info["path_project"], n_cl, dm, alpha_val_DEGs,False)
                    plot_number_pathways_per_gene(DF_sign_pw_per_DEG_gr_i, options["opt_save"], path_info["path_results_subdirectory"], n_cl, path_info["path_colors"], group_name+'_'+add_str)
                
            if options["opt_plot_n_pathways_heatmap"]:
                #heatmaps of number of up-/ down-regulated genes:
                plot_heatmap_pathways_from_DEGs_per_celltype(DF_sign, n_cl, options["opt_save"], path_info["path_results_subdirectory"], variables_groups, variable_n_genes, path_info["path_data"], options["opt_up_and_down_sep"])
                #% genes overlap with geneset
            
            if options["opt_plot_PCA_pathways"]: #up and downregulated together
                if np.shape(DF_sign)[0]>4: #if more than 4 pathways 
                    #this function is very slow:
                    principalDf = get_PCs_for_pathways(DF_sign,variable_gene_names,variables_groups)
                    plot_PCs_pathways(principalDf, n_cl, path_info["path_results_subdirectory"], path_info["path_colors"], 'Celltype', False, variables_groups, options["opt_save"], options["opt_up_and_down_sep"])
                    plot_PCs_pathways(principalDf, n_cl, path_info["path_results_subdirectory"], path_info["path_colors"], variables_groups[0], False, variables_groups, options["opt_save"], options["opt_up_and_down_sep"])
                    
            if options["opt_up_and_down_sep"]:
                DF_all_up, DF_all_down, df_number_pws_not_tested = get_DFs_all_up_down(file_names_all, options["pval_col_to_extract"], alpha_val_DEGs, path_info["path_results_subdirectory"], path_info["path_project"], n_cl, dm, options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"], results_path_proteomics, path_info["path_results_subdirectory"],options["opt_load_dataframe"],options["opt_filter_stats_for_GOs"])
                                        
                DF_all_p_up, DF_all_p_down = get_DFs_as_pivot(DF_all_up, DF_all_down,options["pval_col_to_extract"],options["index_cols"])
                
                if options["opt_plot_pathway_statistics"]:
                    plot_number_of_celltype_specific_and_shared_pathways(pd.concat([DF_all_p_up,DF_all_p_down]),options["p_val_cutoff"],path_info["path_results_subdirectory"], options["opt_save"], add_str, options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],options["opt_filter_stats_for_GOs"])
                #overall stats
                DF_n_pw = get_statistics_pathways(DF_all_p_up, DF_all_p_down, df_number_pws_not_tested, options["p_val_cutoff"], options["pval_col_to_extract"])
            else:
                #proteomics is never seperately tested in this case
                DF_all, df_number_pws_not_tested = get_DFs_all(file_names_all, options["pval_col_to_extract"], alpha_val_DEGs, path_info["path_results_subdirectory"], path_info["path_project"], n_cl, dm, options["opt_proteomics"], results_path_proteomics,path_info["path_results_subdirectory"],options["opt_load_dataframe"],options["opt_filter_stats_for_GOs"])
                                        
                # plot percentage of upregulated DEGs mapping to pathways:
                
                DF_all_p = get_DFs_as_pivot_all(DF_all, options["pval_col_to_extract"],options["index_cols"])
            
                # plot number of cell type specific and shared pathways:
                if options["opt_plot_pathway_statistics"]:
                    plot_number_of_celltype_specific_and_shared_pathways(DF_all_p,options["p_val_cutoff"],path_info["path_results_subdirectory"], options["opt_save"], add_str,options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],options["opt_filter_stats_for_GOs"])
            
                #overall stats
                DF_n_pw = get_statistics_pathways_all(DF_all_p, df_number_pws_not_tested, options["p_val_cutoff"], options["pval_col_to_extract"])
            if options["opt_filter_stats_for_GOs"] and options["opt_up_and_down_sep"]:
                add_str = "_GOs_only"
            else:
                add_str=""
            if options["opt_plot_pathway_statistics"]:
                plot_pathway_statistics(DF_n_pw, options["opt_save"], path_info["path_results_subdirectory"], add_str,n_cl,path_info["path_colors"],options["opt_proteomics"],add_str)
            
            if options["opt_plot_upsets"]:
                if options["opt_up_and_down_sep"]:
                    plot_number_sign_pathways_as_upset(DF_all_p_up, 'up', '', options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],5,options["opt_save"], path_info["path_results_subdirectory"])
                    plot_number_sign_pathways_as_upset(DF_all_p_down, 'down', '', options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],5,options["opt_save"], path_info["path_results_subdirectory"])
                    plot_number_sign_pathways_as_upset(pd.concat([DF_all_p_up,DF_all_p_down]), 'all_separately', '', options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],5,options["opt_save"], path_info["path_results_subdirectory"])
                    plot_number_sign_pathways_as_upset(pd.concat([DF_all_p_up,DF_all_p_down]), 'up_down_together','', options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],5,options["opt_save"], path_info["path_results_subdirectory"])
                else:
                    plot_number_sign_pathways_as_upset(DF_all_p,'all_together','', options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],0,options["opt_save"], path_info["path_results_subdirectory"])
                    
            #stats per group:
            if options["opt_up_and_down_sep"]:
                groups = np.unique(DF_all_p_up.index.get_level_values('group').tolist() + DF_all_p_down.index.get_level_values('group').tolist())
            else:
                groups = DF_all_p.index.get_level_values('group').tolist()
            groups = groups[groups!='gene-ontology']
            #loop over groups:
            for group_name in groups:
                if options["opt_up_and_down_sep"]:
                    #define DF_all_up_group_i,DF_all_down_group_i:
                    DF_all_p_up_group_i = DF_all_p_up.loc[[group_name]].copy()
                    DF_all_p_down_group_i = DF_all_p_down.loc[[group_name]].copy()
                    #DF_all_p_down_group_i = DF_all_p_down[DF_all_p_down.index.get_level_values('group') == group_name].copy()
                    DF_n_pw_group_i = get_statistics_pathways(DF_all_p_up_group_i, DF_all_p_down_group_i, df_number_pws_not_tested, options["p_val_cutoff"], options["pval_col_to_extract"])
                    # plot number of cell type specific and shared pathways:
                    if options["opt_plot_pathway_statistics"]:
                        plot_number_of_celltype_specific_and_shared_pathways(pd.merge(left=DF_all_p_up_group_i,right=DF_all_p_down_group_i,how="right",right_index=True,left_index=True),options["p_val_cutoff"],path_info["path_results_subdirectory"], options["opt_save"], group_name+'_'+add_str,options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"], False)
                    if options["opt_plot_upsets"]:
                        plot_number_sign_pathways_as_upset(DF_all_p_up_group_i, 'up', group_name, options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],1,options["opt_save"], path_info["path_results_subdirectory"])
                        plot_number_sign_pathways_as_upset(DF_all_p_down_group_i, 'down', group_name, options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],1,options["opt_save"], path_info["path_results_subdirectory"])
                        plot_number_sign_pathways_as_upset(pd.concat([DF_all_p_up_group_i,DF_all_p_down_group_i]), 'all_separately',group_name, options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],1,options["opt_save"], path_info["path_results_subdirectory"])
                        plot_number_sign_pathways_as_upset(pd.concat([DF_all_p_up_group_i,DF_all_p_down_group_i]), 'up_down_together',group_name, options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],1,options["opt_save"], path_info["path_results_subdirectory"])
                else:
                    DF_all_p_group_i = DF_all_p.loc[[group_name]].copy()
                    DF_all_group_i = DF_all[DF_all['group']==group_name].copy() 
                    # plot percentage of upregulated DEGs mapping to pathways:
                    plot_percentage_upregulated_DEGs_in_pathways(DF_all_group_i, options["pval_col_to_extract"],options["p_val_cutoff"],path_info["path_results_subdirectory"], options["opt_save"], group_name)
                    # plot number of cell type specific and shared pathways:
                    if options["opt_plot_pathway_statistics"]:
                        plot_number_of_celltype_specific_and_shared_pathways(DF_all_p_group_i,options["p_val_cutoff"],path_info["path_results_subdirectory"], options["opt_save"], group_name+'_'+add_str,options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],False)
                    DF_n_pw_group_i = get_statistics_pathways_all(DF_all_p_group_i, df_number_pws_not_tested, options["p_val_cutoff"], options["pval_col_to_extract"])
                    if options["opt_plot_upsets"]:
                        plot_number_sign_pathways_as_upset(DF_all_p_group_i,'all_together',group_name, options["p_val_cutoff"], options["pval_col_to_extract"],options["index_cols"],1,options["opt_save"], path_info["path_results_subdirectory"])
                if options["opt_plot_pathway_statistics"]:
                    plot_pathway_statistics(DF_n_pw_group_i, options["opt_save"], path_info["path_results_subdirectory"], add_str+'_'+group_name,n_cl,path_info["path_colors"], options["opt_proteomics"],"")
            
            if options["opt_up_and_down_sep"]:
                DF_sorted_filtered = get_DF_sorted_filtered(DF_all_p_up,DF_all_p_down,[],options["p_val_cutoff"],options["index_cols"],options["opt_create_DF_sorted_filtered"], options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"], background_data_proteomics, j, path_info["path_results_subdirectory"], options["opt_up_and_down_sep"])
            else:
                DF_sorted_filtered = get_DF_sorted_filtered(DF_all_p,[],[],options["p_val_cutoff"],options["index_cols"],options["opt_create_DF_sorted_filtered"], options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"], background_data_proteomics, j, path_info["path_results_subdirectory"], options["opt_up_and_down_sep"])
    return DF_sorted_filtered

# plot percentage of upregulated DEGs mapping to pathways:
def plot_percentage_upregulated_DEGs_in_pathways(DF_all,pval_col_to_extract,p_val_cutoff, path_results_long, opt_save, add_str):
    plt = loadPltSettings(8,1)
    plt.ticklabel_format(style='plain', axis='x',useOffset=False)
    fig, ax = plt.subplots(figsize=(9,6))
    #only significant pathways
    #TypeError: '<=' not supported between instances of 'str' and 'float'
    DF_all[pval_col_to_extract] = DF_all[pval_col_to_extract].str.replace(',','.').astype(float)
    DF_sign = DF_all[DF_all[pval_col_to_extract]<=p_val_cutoff].copy()
    DF_sign['Percentage upregulated genes'] = 100*DF_sign['number upregulated DEGs']/(DF_sign['number upregulated DEGs']+DF_sign['number downregulated DEGs'])
    if len(DF_sign)>0:
        sns.stripplot(data=DF_sign,y='celltype_name', x='Percentage upregulated genes',ax=ax,color='white',edgecolor='black',size=3)
        
        y1,y2 = ax.get_ylim()
        yrange=math.ceil(np.abs(y2-y1))
        y_start = min(y1,y2)
        x1,x2 = ax.get_xlim()
        xrange=100
        ax.add_patch(Rectangle((0, y_start),xrange/3,yrange,facecolor = 'teal',fill=True))
        ax.add_patch(Rectangle((xrange/3, y_start),xrange/3,yrange,facecolor = 'grey',fill=True))
        ax.add_patch(Rectangle((2*xrange/3, y_start),xrange/3,yrange,facecolor = 'darkred',fill=True))
        ax.ticklabel_format(style='plain',axis='x')
        ax.set_xlim([0,100])
        
        if opt_save:
            path=path_results_long+'figures/stats/'
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            fig.savefig(path + 'percentage_upregulated_DEGs_in_pathways_per_celltype'+add_str.replace(':','_')+'.pdf',bbox_inches='tight')
            plt.close(fig)
            plt.clf()   

def plot_number_sign_pathways_as_upset(DF, r_str, group_name, p_val_cutoff, pval_col_to_extract,index_cols,min_ss,opt_save, path_results):
    #TO DO: short ct names
    params = {'font.size':10, 'axes.labelsize':10, 'xtick.labelsize':10, 'ytick.labelsize':10}
    if r_str == 'up':
        color='darkred'
    elif r_str =='down':
        color='teal'
    else:
        color='grey'
    DF = get_short_CT_names_as_columns(DF)
    CTs = DF.columns.tolist()
    #check if columns have correct data type to create bool df:
    if not np.all(DF.dtypes==float):
        print("To DO: convert columns!")
    if r_str == 'up_down_together':
        #create DF containing bools if p-vals are lower than cutoff:
        DF_bool = DF.lt(p_val_cutoff)
        #all CTs:
        #To do change as they do end with " up" and " down" now instead of "_up" and "_down"
        CTs_raw = np.unique([item.split(' up', 1)[0] for item in CTs if item.endswith('up')] + [item.split(' down', 1)[0] for item in CTs if item.endswith('down')])
        #sum up up- and down-reg:
        for ct in CTs_raw.tolist():
            if ct+' up' in DF.columns and ct+' down' in DF.columns:
                DF_bool[ct] = DF_bool[ct+' up'] | DF_bool[ct+' down']
        DF_bool.drop(columns=CTs,inplace = True)
        CTs=CTs_raw
    else:
        #create DF containing bools if p-vals are lower than cutoff:
        DF_bool = DF.lt(p_val_cutoff)
    # # remove all that are False for all CTs and proteomics:
    DF_bool_red = DF_bool[DF_bool.sum(axis=1)>0].copy()
    
    #min subset size in upset plot
    #at least one CT needs to have more than min_ss
    # at least one pathway needs to be implicated in more than one cell type
    if np.max(DF_bool_red.sum(axis=0))>min_ss and np.any(DF_bool_red.sum(axis=1)>1):
        #rename columns (add spaces)
        column_names = DF_bool_red.columns.tolist()
        for col in column_names:
            eval("DF_bool_red.rename(columns={\'"+col+"\':\'               "+col+"\'},inplace=True)")
        # update CTs
        CTs = DF_bool_red.columns.tolist()
        if len(CTs)>np.max([2,min_ss]):
            with plt.rc_context(params):
                plot(from_indicators(CTs,data=DF_bool_red), show_counts=True, sort_by = 'cardinality', sort_categories_by='cardinality', min_subset_size=min_ss, facecolor=color)# 
            if opt_save:
                if group_name=='':
                    path=path_results+'figures/upset/'
                else:
                    path=path_results+'figures/upset/'+group_name.replace(':','_')+'/'
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                plt.savefig(path+'upset_'+r_str+'_significant_pathways_'+str(p_val_cutoff).replace('.','_')+'.pdf', bbox_inches='tight')
            plt.close()
            plt.clf() 

def plot_pathway_statistics(DF_n_pw,opt_save, path_results_long, group_name, n_clusters,path_filtered_data,opt_proteomics,add_str):

    CT_ordered_m = get_CT_order(n_clusters,path_filtered_data,True)
    #CT_ordered_m = [ct.replace(' ','_').replace('-','_') for ct in CT_ordered]

    if opt_proteomics:
        CT_ordered_m.append('proteomics')
        palette_dict = {'up':'darkred','down':'teal','deregulated': 'navy'}
        hue_vec = ['up','down','deregulated']
    else:
        palette_dict = {'up':'darkred','down':'teal'}
        hue_vec = ['up','down']
    
    plt = loadPltSettings(8,1)
    #to do: xtick labels CTs not classes
    #legend: show up-reg down-reg de-reg, not which CT
    #to do: distinguish between proteomics up and down (1 subplot) and proteomics all (2 subplots with shared axes)
    # plts.use('ggplot')
    fig, ax = plt.subplots(3,1,figsize=(5,9), tight_layout=True, sharex=True)
    fig.suptitle('Pathway stats '+group_name)
    sns.barplot(data=DF_n_pw,x='Celltype_short',y='significant',hue='regulation',hue_order=hue_vec,palette=palette_dict,ax=ax[0],order=CT_ordered_m, edgecolor='white')
    ax[0].set_xticklabels(ax[0].get_xticklabels() , rotation=90)
    ax[0].legend(loc='upper right')
    ax[0].ticklabel_format(style='plain', axis='y')
    ax[0].set(ylabel='Significant pathways')
    sns.barplot(data=DF_n_pw,x='Celltype_short',y='Percentage significant',hue='regulation',hue_order=hue_vec,palette=palette_dict,ax=ax[1],order=CT_ordered_m, edgecolor='white')
    ax[1].set_xticklabels(ax[1].get_xticklabels() , rotation=90)
    ax[1].set_ylim([0,100])
    ax[1].legend([],[], frameon=False)
    ax[1].ticklabel_format(style='plain', axis='y')
    ax[1].set(ylabel='% significant pathways')
    #plot number of pathways tested:
    sns.barplot(data=DF_n_pw,x='Celltype_short',y='Percentage tested',hue='regulation',hue_order=hue_vec,palette=palette_dict,ax=ax[2],order=CT_ordered_m, edgecolor='white')
    ax[2].set_xticklabels(ax[2].get_xticklabels() , rotation=90)
    ax[2].set_ylim([0,100])
    ax[2].legend([],[], frameon=False)
    ax[2].ticklabel_format(style='plain', axis='y')
    ax[2].set(xlabel='Celltype', ylabel='% pathways tested')
    if opt_save:
        path=path_results_long+'figures/stats/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        fig.savefig(path + 'sign_pathways_'+group_name.replace(':','_')+add_str+'.pdf',bbox_inches='tight')
        plt.close(fig)
        plt.clf()   
   
def plot_PCs_pathways(principalDF, n_cl, path_results, path_colors, opt_color_code, opt_with_labels,variables_groups,opt_save, opt_transcriptomics_up_and_down_sep):

    plt = loadPltSettings(14,3)
    if opt_color_code =='Celltype':
        #define color palette:

        if "Patient_clustering/" in path_results:
            if opt_transcriptomics_up_and_down_sep:
                principalDF = principalDF.iloc[~principalDF.index.get_level_values('Celltype').str.endswith('__up_or_down_sig')]
            else:
                principalDF = principalDF.iloc[principalDF.index.get_level_values('Celltype').str.endswith('__up_or_down_sig')]
        _,idx = np.unique(principalDF.index.get_level_values("Celltype"),return_index=True)
        CTs = principalDF.index.get_level_values("Celltype")[np.sort(idx)].tolist() # preserved order!
        if opt_transcriptomics_up_and_down_sep:
            if "Patient_clustering/" in path_results:
                CTs = [ct.replace('__up_sig','').replace('__down_sig','') for ct in CTs]
            else:
                CTs = [ct.replace('_up_sig','').replace('_down_sig','').replace('_',' ').replace(' 2 3 ', ' 2-3 ').replace(' 3 4 ', ' 3-4 ').replace(' 3 6 ', ' 3-6 ').replace(' 5 6 ', ' 5-6 ') for ct in CTs]
        else:
            if "Patient_clustering/" in path_results:
                CTs = [ct.replace('__up_or_down_sig','') for ct in CTs]
            else:
                CTs = [ct.replace('_sig','').replace('_',' ').replace(' 2 3 ', ' 2-3 ').replace(' 3 4 ', ' 3-4 ').replace(' 3 6 ', ' 3-6 ').replace(' 5 6 ', ' 5-6 ') for ct in CTs]
        #search for color in sheet:
        os.chdir(path_colors)
        df_colors = pd.read_excel('Cell_type_colors.xlsx',engine="openpyxl")
        palette=[]
        if n_cl == 3:
            column_ct_i = 'Cell class ('+str(n_cl)+')'
            column_color_i = 'Cell class ('+str(n_cl)+') color'
        else:
            column_ct_i = 'Cell type ('+str(n_cl)+')'
            column_color_i = 'Cell type ('+str(n_cl)+') color'
        CTs_color_sheet = np.unique(df_colors[column_ct_i].tolist())
        for ct_i in CTs:
            if ct_i in CTs_color_sheet:
                palette = palette + np.unique(df_colors[df_colors[column_ct_i]==ct_i][column_color_i]).tolist()

        fig, ax = plt.subplots()
        sns.scatterplot(x=principalDF['principal component 1'] , y=principalDF['principal component 2'] , hue=principalDF.index.get_level_values("Celltype"), palette = palette, legend='full', ax=ax).set(title='pathways of deregulated genes')
    else:
        fig, ax = plt.subplots()
        sns.scatterplot(x=principalDF['principal component 1'] , y=principalDF['principal component 2'] , hue=principalDF.index.get_level_values(variables_groups[0]), legend='full', ax=ax).set(title='pathways of deregulated genes')
    if opt_with_labels:
        if np.shape(principalDF)[0]<40:
            texts = []
            if np.shape(principalDF)[0]>20:
                principalDF['go_pathway']=principalDF.index.get_level_values(variables_groups[1])
            else:
                principalDF['go_pathway']=principalDF.index.get_level_values(variables_groups[0])+'__'+principalDF.index.get_level_values(variables_groups[1])
            for r, s in principalDF[['principal component 1','principal component 2', 'go_pathway']].iterrows():
                texts.append(plt.text(s['principal component 1'],s['principal component 2'],s['go_pathway'],fontSize=3))
                #make sure labels are non-overlapping:
                adjust_text(texts, only_move={'points':'y', 'texts':'y'})
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if opt_save:
        path=path_results+'figures/PCA/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'pathways_GSA_deregulated_colored_by_'+opt_color_code+'.pdf',bbox_inches='tight')
        plt.close()
        plt.clf()

### plot functions
def plot_signed_DEG_corr_clustermaps(df_DESeq2_all_genes_all_CTs,path_results):
    df_DESeq2_all_genes_all_CTs = get_short_CT_names(df_DESeq2_all_genes_all_CTs,ct_column="celltype")
    df_DESeq2_all_genes_all_CTs["signed_padj"] = df_DESeq2_all_genes_all_CTs["padj"]
    df_DESeq2_all_genes_all_CTs.loc[df_DESeq2_all_genes_all_CTs["log2FoldChange"]<0,"signed_padj"] = df_DESeq2_all_genes_all_CTs[df_DESeq2_all_genes_all_CTs["log2FoldChange"]<0]["padj"]*(-1)
    pivot_opt = ["signed_padj","log2FoldChange"]
    for ps in pivot_opt:
        df_tmp = df_DESeq2_all_genes_all_CTs[[ps,"celltype_short"]].copy()
        df_pivot = pd.pivot_table(data=df_tmp,values=ps,columns="celltype_short",index="Gene_short")
        corr_methods = ["pearson","spearman"]
        for m in corr_methods:
            df_corr = df_pivot.corr(method=m)
            sns.clustermap(data=df_corr,cmap="BrBG",vmin=-1, vmax=1)
            path=path_results+'figures/correlation_between_CTs/'
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            plt.savefig(path+'Genes_'+m+'_correlation_between_cts_based_on_'+ps+'.pdf',bbox_inches='tight')
            plt.close()
            plt.clf()

def plot_heatmap_DEG_result_for_genes_in_genelist(DF_DEG,genelist,genelist_name,q_values,opt_save,path_results,opt_genes_to_show,opt_par,genelist_version,mode,opt_add_proteomics,path_results_proteomics, path_results_proteomics_hgnc, hgnc_file, proteomics_file):
    if genelist_name in ['asd.wes','scz.wes']:
        string1 = 'q-value '
    else:
        string1 = 'p-value '
    qval_cutoff=0.001
    if 'wes' in genelist_name:
        n_round = 4
        #check if qvalues empty:
        if len(q_values)>0:
            qval_cutoff = 0.1
    else:
        n_round=7
        
    if len(q_values)>0:
        #filter DF for genes in genelist
        #make sure order of genes is the same
        #best way: create another dataframe and merge them 
        DF = get_genelist_DEG_DF(DF_DEG, genelist, genelist_name, q_values, string1)
        #reduce to q-values below qval_cutoff:
        DF_red_all = DF[DF[string1+genelist_name]<qval_cutoff].copy()
        DF_red_all.sort_values(by=string1+genelist_name,inplace=True)
        _,idx =np.unique(DF_red_all['Gene'], return_index=True)
        opt_q_val_in_gl=True
    else:
        #print('what happens now?')
        DF_red_all = get_genelist_DEG_DF(DF_DEG, genelist, genelist_name, q_values,string1)
        idx=range(0,len(DF_red_all))
        opt_q_val_in_gl=False
    
    #add line for proteomics in graphic
    if opt_add_proteomics:
        DF_p = get_DEP_results(path_results_proteomics,proteomics_file,1,'','', 1 ,path_results_proteomics_hgnc+hgnc_file)
        #add toDF_red_all as separate "celltype"
        if 'celltype_short' in DF_red_all.columns:
            ct_col = 'celltype_short'
        elif 'celltype_name' in DF_red_all.columns:
            ct_col = 'celltype_name'
        else:
            ct_col = 'celltype'
        DF_p[ct_col] = "Proteomics"
        DF_p=DF_p.rename(columns={"Gene":"Gene_short","log2FoldChange_DEP":"log2FoldChange","qvalue_DEP":"padj"})
        DF_p.drop(columns=[c for c in DF_p.columns if c not in ["Gene_short","log2FoldChange","padj",ct_col]],inplace=True)#,"signed_padj","Color_sign"]])
        #only keep genes of interest (in current genelist):
        genes_of_interest = np.unique(DF.index.tolist()).tolist()
        DF_red_all = DF_red_all.reset_index()
        #merge data frames:
        DF_red_all = pd.concat([DF_red_all,DF_p])
        DF_red_all = DF_red_all.set_index("Gene_short")
        gene_sel = [g for g in np.unique(DF_red_all.index.tolist()).tolist() if g in genes_of_interest]
        DF_red_all = DF_red_all.loc[gene_sel]

    if opt_genes_to_show=='top_genes':
        #keep only top genes:
        genes_shortlisted = DF_red_all['Gene'][np.sort(idx)].tolist()
        genes_shortlisted_unique = get_unique_values_in_list_while_preserving_order(genes_shortlisted)
        number_of_genes_total = len(genes_shortlisted_unique)
        if number_of_genes_total>opt_par:
            rounds = round(number_of_genes_total/opt_par)
            genes_of_interest = genes_shortlisted_unique[0:opt_par]
        else:
            rounds = 1
            genes_of_interest = genes_shortlisted_unique
        for r_id in range(0,rounds):
            if r_id>0 and opt_genes_to_show=='top_genes':
                genes_of_interest = genes_shortlisted_unique[opt_par*(r_id-1):opt_par*r_id]
            DF_red = DF_red_all[DF_red_all['Gene'].isin(genes_of_interest)]
            heatmap_data, _, number_of_genes_plotted = get_ordered_heatmap_data_genelist_DEG_overlap_and_number_of_genes(DF_red, genelist_name, n_round, opt_par, opt_genes_to_show,opt_q_val_in_gl,mode,opt_add_proteomics)
            if number_of_genes_plotted > 0:
                # _,idx = np.unique(DF_red['Gene (q_val)'],return_index=True)
                # idx_sorted = DF_red['Gene (q_val)'][np.sort(idx)]
                # heatmap_data = heatmap_data.reindex(idx_sorted)
                #color_def = np.array(fc.dec_to_rgb(colors.to_rgb('white'))).tolist(), np.arrayfc.dec_to_rgb(colors.to_rgb('teal')).tolist(),np.array(fc.dec_to_rgb(colors.to_rgb('darkred')).tolist(),np.array(fc.dec_to_rgb(colors.to_rgb('white')).tolist()
                plot_heatmap_genelist(heatmap_data, r_id, mode, DF_red_all, genelist_name,opt_genes_to_show,opt_save, path_results,opt_add_proteomics)
            
    elif opt_genes_to_show=='genes_with_strong_signal':
        genes_of_interest = DF_red_all['Gene'][np.sort(idx)].tolist()            
        DF_red = DF_red_all[DF_red_all['Gene'].isin(genes_of_interest)]
        if np.shape(DF_red)[0]>0:
            print(genelist_name)
            heatmap_data, number_of_genes_total, number_of_genes_plotted = get_ordered_heatmap_data_genelist_DEG_overlap_and_number_of_genes(DF_red, genelist_name, n_round, opt_par, opt_genes_to_show,opt_q_val_in_gl,mode,opt_add_proteomics)
            # _,idx = np.unique(DF_red['Gene (q_val)'],return_index=True)
            # idx_sorted = DF_red['Gene (q_val)'][np.sort(idx)]
            # heatmap_data = heatmap_data.reindex(idx_sorted)
            #color_def = np.array(fc.dec_to_rgb(colors.to_rgb('white'))).tolist(), np.arrayfc.dec_to_rgb(colors.to_rgb('teal')).tolist(),np.array(fc.dec_to_rgb(colors.to_rgb('darkred')).tolist(),np.array(fc.dec_to_rgb(colors.to_rgb('white')).tolist()
            if number_of_genes_plotted > 0:
                rounds = int(np.ceil(np.shape(heatmap_data)[0]/50))
                for r_id in range(0,rounds):
                    heatmap_data_sel = heatmap_data.iloc[50*r_id:50*(r_id+1), ]
                    plot_heatmap_genelist(heatmap_data_sel, r_id, mode, DF_red_all, genelist_name,opt_genes_to_show,opt_save, path_results,opt_add_proteomics)
        else:
            number_of_genes_total = 0
            number_of_genes_plotted = 0
    return (number_of_genes_total, number_of_genes_plotted)

def plot_heatmap_genelist(heatmap_data, r_id, mode, DF_red_all, genelist_name,opt_genes_to_show,opt_save, path_results):
    
    if mode == 'p-value':
        new_cmap = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('teal'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('darkred'))).tolist())#([64,57,144],[112,198,162],[230,241,146],[253,219,127],[244,109,69],[169,23,69])
    elif mode=='log2FC':
        new_cmap = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('teal'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('darkred'))).tolist())
    new_cmap.set_bad("silver") 
    
    
    if mode == 'p-value':
        ax = sns.heatmap(heatmap_data.transpose(), yticklabels=1, xticklabels=1, vmin = -1, vmax = 1, cmap = new_cmap, square=True, cbar_kws={'ticks' : [-1,-0.9,-0.8,0,0.8,0.9,1],'label': 'adjusted ' +mode}, mask=heatmap_data.transpose().isnull())
    else:
        max_LFC_abs = np.max([np.abs(DF_red_all['log2FoldChange'].min()), np.abs(DF_red_all['log2FoldChange'].max())])
        ax = sns.heatmap(heatmap_data.transpose(), yticklabels=1, xticklabels=1, vmin = (-1)*max_LFC_abs, vmax = max_LFC_abs, cmap = new_cmap, square=True, mask=heatmap_data.transpose().isnull(),cbar_kws={'ticks' : [(-1)*max_LFC_abs,0,max_LFC_abs],'label': mode})
    ax.figure.axes[-1].yaxis.label.set_size(12)
    ax.axhline(sum(heatmap_data.columns.str.startswith('Exc')),c="grey",linewidth=1)
    ax.axhline(sum(heatmap_data.columns.str.startswith('Exc'))+sum(heatmap_data.columns.str.startswith('Inh')),c="grey",linewidth=1)
    
    
    if np.shape(heatmap_data)[0]>30:
        plt.yticks(fontsize=7)
        plt.xticks(fontsize=7)
    else:
        plt.yticks(fontsize=10)
        plt.xticks(fontsize=10)
            
    if mode == 'p-value':
        ax.collections[0].colorbar.set_ticklabels(['0 (down)','0.1','0.2','1','0.2','0.1','0 (up)'])
    ax.set_title(genelist_name, size=12)
    #Below 3 lines remove default labels
    ax.set_xlabel('Gene in ' + genelist_name, fontsize=12)
    ax.set_ylabel('Celltype', fontsize=12)

    if opt_save:
        path=path_results+'figures/comparison_with_genelists/'+opt_genes_to_show+'/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'heatmap_'+genelist_name.replace('.','_')+'_'+str(r_id)+'_'+mode+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf() 

def plot_percentage_genes_in_list_and_DEG_for_range_of_cutoffs(path_project, DF_DEG, opt_par, qval_cutoff_range, path_results_long, opt_save): 
    P_genes, genelist_names = get_percentage_genes_in_list_and_DEG_for_range_of_cutoffs(path_project, DF_DEG, opt_par, qval_cutoff_range)
    #define colors:
    cmap = plt.get_cmap('autumn')
    colors = cmap(np.linspace(0, 1, len(qval_cutoff_range)))
    #sort genelist_names and P_genes according to lowest qval_cutoff:
    vec_unsorted = P_genes[:,0]
    idx_sorted = [x for x, y in sorted(enumerate(vec_unsorted), key=lambda x: x[1],reverse=True)]
    P_genes_sorted = P_genes[idx_sorted,:]
    genelist_names_sorted = [genelist_names[idx] for idx in idx_sorted]   
    #remove list that appears twice:
    idx_to_plot = [idx for idx, element in enumerate(genelist_names_sorted) if element!='list_scz.wes']
    idx_remove = [idx for idx, element in enumerate(genelist_names_sorted) if element=='list_scz.wes']
    del genelist_names_sorted[idx_remove[0]]
    
    #plot stacked barplot
    w=0.6 # width of stacked chart
    handles = []
    for j in range(np.shape(P_genes_sorted)[1]-1,-1,-1):
        if j==0:
            add_str = '0 <='
            c = 'darkred'
        else:
            add_str = str(qval_cutoff_range[j-1]) + ' < '
            c = colors[j]
        h = plt.bar(genelist_names_sorted,P_genes_sorted[idx_to_plot,j],w,label = add_str + 'p_adj <= '+ str(qval_cutoff_range[j]),color = c,edgecolor="white")
        handles.append(h)
    plt.ticklabel_format(axis="y", style='plain')
    plt.xticks(rotation='vertical')
    plt.ylabel('Percentage of genes')
    plt.legend(bbox_to_anchor=(1.005, 1.02), loc='upper left')
    if opt_save:
        path=path_results_long+'figures/comparison_with_genelists/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'Percentage_genes_overlapping_with_genelists.pdf', bbox_inches='tight')
        plt.close()
        plt.clf()

def plot_log2FC_per_CT(df_DESeq2_all_genes,alpha_val,n_clusters,path_filtered_data,path_results_long,opt_save):
    #add short CT labels
    nan = np.nan
    #get short ct names
    df_DESeq2_all_genes['celltype'] = df_DESeq2_all_genes['celltype'].str.replace('_',' ')
    df_DESeq2_all_genes = get_short_CT_names(df_DESeq2_all_genes,"celltype")

    #make sure CT are correctly ordered
    if n_clusters!=3:
        CTs_sorted = get_CT_order(n_clusters,path_filtered_data,True)
    else:
        CTs_sorted = np.unique(df_DESeq2_all_genes['celltype']).tolist()
    #
        #otherwise names don't match
    df_DESeq2_all_genes=df_DESeq2_all_genes.sort_values('celltype_short', key=make_sorter(CTs_sorted))

    #ax = sns.boxplot(y=df_DESeq2_all_genes['log2FoldChange'],x=df_DESeq2_all_genes['celltype'],color='lightgray',linewidth=1,flierprops = dict(markerfacecolor = '1', markersize = 3,marker='o'))
    plt = loadPltSettings(12,2)
    fig, ax = plt.subplots(figsize=(9,6))
    sns.violinplot(data = df_DESeq2_all_genes, y='log2FoldChange',x='celltype_short',color='lightgray',linewidth=1,ax = ax)
    plt.setp(ax.get_xticklabels(), rotation=90)
    CTs = np.unique(df_DESeq2_all_genes['celltype_short']).tolist()
   
    reg_col_str = ['darkred','teal']
    for rcs in reg_col_str:
        DF_sel = df_DESeq2_all_genes[df_DESeq2_all_genes.Color_sign==rcs].copy()
        for ct in CTs:
            #check which CTs are missing and add dummy variable:
            if np.sum(DF_sel['celltype_short']==ct)==0:
                DF_sel = pd.concat([DF_sel, pd.DataFrame(data={'Gene':nan, 'baseMean':nan, 'log2FoldChange':nan, 'lfcSE':nan, 'stat':nan, 'pvalue':nan, 'padj':nan, 'negative_log10_padj':nan, 'Color_sign':nan, 'celltype_short':ct},index=[0])],ignore_index=True)
        DF_sel=DF_sel.sort_values('celltype_short', key=make_sorter(CTs_sorted))
        sns.stripplot(data=DF_sel, x="celltype_short", y="log2FoldChange", ax=ax, color=rcs,size=2, zorder=1, dodge=True)
    ax = remove_frame_keep_axes(ax)
    ax.set_ylim([-2,2])
    
    # plt.xlim(-0.5, len(np.unique(df_DESeq2_all_genes['celltype']))-0.5)
    fig = ax.get_figure()
    if opt_save:
        fig.savefig(path_results_long + 'figures/Log2FC_per_CT'+str(alpha_val).replace('.','_')+'.pdf',bbox_inches='tight')
        plt.close(fig)
        plt.clf()   

#form https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
def make_sorter(l):
    """
    Create a dict from the list to map to 0..len(l)
    Returns a mapper to map a series to this custom sort order
    """
    sort_order = {k:v for k,v in zip(l, range(len(l)))}
    return lambda s: s.map(lambda x: sort_order[x])
   
def plot_cells_per_donor_for_each_cell_type_cluster(df,df_donor_disease,df_aggregation_summary,cutoff_p,opt_save,path_results,opt_mode,opt_celltype_groups,n_clusters,cl_str):
    #plot number of cells/ percentages of cells per donor for each cell type cluster
    G=['CTRL','SCZ']
    colours = {"CTRL": '#487BAF', "SCZ": '#FA9835'}
    plt = loadPltSettings(10,20)
    if n_clusters == 3:
        n_rows = 1
    elif  n_clusters == 15:
        n_rows = 7
    else:
        n_rows= 18
    fig, axes = plt.subplots(n_rows,  3, figsize=(45,55)) 
    fig.subplots_adjust(hspace = .3, wspace = 0.1)  
    counter_N =-1
    counter_E=-1
    counter_I=-1

    CT_cols = [ct for ct in df.columns.tolist() if ct not in ['Donor','Disease']]
    
    for idx,col in enumerate(CT_cols):
        if col.startswith('NonNeuronal') | col.startswith('Astr') | col.startswith('Olig') | col.startswith('Micro'):
            col_idx = 0
            counter_N = counter_N +1
            row_idx = counter_N
        elif col.startswith('Excitatory'):
            col_idx = 1
            counter_E = counter_E+1
            row_idx=counter_E
        elif col.startswith('Inhibitory'):
            col_idx = 2   
            counter_I = counter_I+1
            row_idx=counter_I
        else:
            continue
        print("counter_N: "+ str(counter_N))
        print("counter_E: " + str(counter_E))
        print("counter_I: " + str(counter_I) )
        if n_clusters==3:
            AX = axes[col_idx]
        else:
            AX = axes[row_idx,col_idx]
        if n_clusters==3:
            df.plot.bar(x='Donor',y=col,color=df['Disease'].replace(colours),title=col,ax=AX)
            AX.set_ylabel(opt_mode+' of cells')
            AX.legend(G)
            AX.plot([0,np.shape(df_donor_disease)[0]],[100/np.shape(df_donor_disease)[0],100/np.shape(df_donor_disease)[0]],'--k')
        y_val_line=[]
        if not df_aggregation_summary.empty:
            y_val_line = df_aggregation_summary[df_aggregation_summary['Cell type cluster']==col]['n_cells_selected'].tolist()[0]
        if cutoff_p:
            y_val_line = cutoff_p
        if y_val_line:
            if n_clusters==3:
                AX.plot([0,np.shape(df_donor_disease)[0]],[y_val_line,y_val_line],'-r')
    if n_clusters!=3:
        if counter_N+1<n_rows:
            for r_id_del in range(counter_N+1,n_rows):
                fig.delaxes(axes[r_id_del,0])
        if counter_E+1<n_rows:
            for r_id_del in range(counter_E+1,n_rows):
                fig.delaxes(axes[r_id_del,1])
        if counter_I+1<n_rows:
                for r_id_del in range(counter_I+1,n_rows):
                    fig.delaxes(axes[r_id_del,2])
    if opt_save:
        fig.savefig(path_results+opt_mode+'_cells_per_Sample_and_CT_cluster_'+opt_celltype_groups+cl_str+'.pdf', bbox_inches='tight')
    plt.close(fig)
    plt.clf()    


def plot_number_of_deregulated_genes(DF_DEG_results,shrinkage_method,opt_save,path_results,opt_DEG_method,add_figure_str,option,n_cl, path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_data,gene_matrix_file,n_clusters,opt_angel_plot):
    # to do: distinguish in protein coding and non-protein coding
    if option == 'vs_number_cells':
        plt = loadPltSettings(15,35)
    else:
        if n_cl >20:
            plt = loadPltSettings(12,35)
        else:
            plt = loadPltSettings(20,35)
    #apply alpha threshold:
    DF_DEG_results = DF_DEG_results[DF_DEG_results["padj"]<=alpha_val]    
    
    #fig, axes = plt.subplots(3,  1, figsize=(20,10), sharex=True) 
    #fig.subplots_adjust(hspace = .3, wspace = 0.1)  
    if np.shape(DF_DEG_results)[0]>0:
        if opt_DEG_method == 'Wilcoxon':
            DF_number = get_number_of_deregulated_genes(DF_DEG_results, 'logfoldchanges', 'Cell type', 'names', n_cl, path_filtered_data)
        elif opt_DEG_method=='DESeq2':
            if 'shrinkage' in DF_DEG_results.columns:
                if ('Unnamed: 0' in DF_DEG_results.columns):
                    DF_DEG_results[["genes","ensgid"]] = DF_DEG_results["Unnamed: 0"].str.split(pat='_',n=1,expand=True)
                    DF_DEG_results.drop(columns=['Unnamed: 0'],inplace=True)
                DF_number = get_number_of_deregulated_genes(DF_DEG_results[DF_DEG_results['shrinkage']==shrinkage_method].copy(), 'log2FoldChange', 'celltype', 'genes', n_cl, path_filtered_data)
                if option != 'vs_number_cells':
                    DF_number = get_short_CT_names(DF_number.reset_index(),'celltype')
                    DF_number.set_index("celltype_short",inplace = True)
                    DF_number.index.name = 'Cell type'
                    #DF_number.drop(columns = ["celltype"],inplace=True)
            else:
                # if 'celltype_name' in DF_DEG_results.columns:
                #     ct_col = 'celltype_name'
                # else:
                if "celltype" in DF_DEG_results.columns:
                    ct_col = "celltype" #'celltype_names'
                elif "celltype_names" in DF_DEG_results.columns:
                    ct_col = "celltype_names"
                else:
                    ct_col = "no clue"
                DF_number = get_number_of_deregulated_genes(DF_DEG_results, 'log2FoldChange', ct_col, 'Gene', n_cl, path_filtered_data)
                if option != 'vs_number_cells':
                    DF_number = get_short_CT_names(DF_number.reset_index(),ct_col)
                    DF_number.set_index(ct_col+"_short",inplace = True)
                    DF_number.index.name = ct_col
                    DF_DEG_results[["genes","ensgid"]] = DF_DEG_results["Gene"].str.split(pat='_',n=1,expand=True)
                    DF_DEG_results_detailed = get_gene_type_for_deregulated_genes(DF_DEG_results, 'log2FoldChange', gene_matrix_file, ct_col)
                    DF_DEG_results_detailed = get_short_CT_names(DF_DEG_results_detailed.reset_index(),ct_col)
                    DF_DEG_results_detailed.set_index(ct_col+"_short",inplace = True)
                    #plot_donut_chart_DEGs(DF_DEG_results_detailed,opt_save, path_results,opt_DEG_method,shrinkage_method, option,opt_celltype_groups,add_figure_str, alpha_val,n_clusters,path_filtered_data)
                    plot_detailed_stacked_bar_DEGs(DF_DEG_results_detailed,opt_save, path_results,opt_DEG_method,shrinkage_method, option,opt_celltype_groups,add_figure_str, alpha_val,n_clusters,path_filtered_data)

        else: 
            print("DEG_method "+opt_DEG_method+" is not specified!")
        DF_number.index = DF_number.index.str.replace('_',' ')
        if option == 'vs_celltypes':
            #remove underscore in Cell type names in DF_number
            if opt_angel_plot and n_cl==16:
                #get enrichment data:
                DF_e = pd.read_csv(path_results+"enrichment_per_CT_dataframe_p3.LDSC_50kb_autosome.tsv",sep="\t")
                DF_e = get_short_CT_names(DF_e,"VARIABLE")
                DF_e["-log10(p-value)"] = np.log10(DF_e["P"])
                DF_e.set_index("VARIABLE_short",inplace=True)
                DF_both = pd.merge(left = DF_e,right=DF_number,left_index=True,right_index=True).sort_values(by="total",ascending=False)
                fig, axes = plt.subplots(1,  2, figsize=(12,5), sharey=True, gridspec_kw={'wspace': 0.0, "width_ratios":[3,2]}) 
                DF_both[['down-regulated','up-regulated']].plot(kind='barh',align="center",stacked=True,xlabel='Number of genes',color={"down-regulated":"teal","up-regulated":'darkred'},edgecolor="white",ax=axes[0])
                axes[0].invert_yaxis()
                axes[0].invert_xaxis()
                axes[0].ticklabel_format(style='plain', axis='x') 
                axes[0] = remove_frame_keep_axes(axes[0])
                DF_both['-log10(p-value)'].plot(kind='barh', align="center",stacked=False,xlabel='-log10(p-value)',color=DF_both["ct_color"],edgecolor="white",ax=axes[1])
                axes[1].axvline(x=np.log10(0.05/np.shape(DF_both)[0]),ls="--",color="grey")
                axes[1].invert_yaxis()
                axes[1].invert_xaxis()
                axes[1] = remove_frame_keep_axes(axes[1])
            else: 
                axes=DF_number[['down-regulated','up-regulated']].plot(kind='barh',stacked=True,ylabel='Number of genes',color={"down-regulated":"teal","up-regulated":'darkred'},edgecolor="white")
                axes.ticklabel_format(style='plain', axis='x') 
                plt.gca().invert_yaxis()
                axes = remove_frame_keep_axes(axes)
                #DF_DEG_results_up = DF_DEG_results[DF_DEG_results['logfoldchanges']>0]
                #DF_DEG_results_down = DF_DEG_results[DF_DEG_results['logfoldchanges']<0]
                #DF_DEG_results.set_index(['Cell type','names']).count(level='Cell type')['scores'].sort_values().plot(kind='bar',title='Deregulated genes',ylabel='Genes',color='grey',ax=axes[0])
                #DF_DEG_results_up.set_index(['Cell type','names']).count(level='Cell type')['scores'].sort_values().plot(kind='bar',title='Upregulated genes',ylabel='Genes',color='teal',ax=axes[1])
                #DF_DEG_results_down.set_index(['Cell type','names']).count(level='Cell type')['scores'].sort_values().plot(kind='bar',title='Downregulated genes',ylabel='Genes',color='sandybrown',ax=axes[2])
        elif option == 'vs_number_cells':
            filename_DF_number = "DF_number_cells_"+str(n_cl)+"_number_DEGs_"+str(alpha_val).replace('.','')+".csv"
            if os.path.isfile(path_data+filename_DF_number):
                DF_number = pd.read_csv(path_data+filename_DF_number,index_col=0)
            else:
                os.chdir(path_filtered_data)
                loom_filename = "Samples_"+opt_celltype_groups[:5]+add_str+add_str_pagoda+"_TH_and_D_adj__filtered_and_CT_annotated_and_CT_clustered.loom"
                n_cells=[]
                with loompy.connect(loom_filename) as D:
                    #add number of cells for each cluster
                    if n_cl == 3:   
                        CT_list = DF_number.index.tolist()#['Excitatory','Inhibitory','NonNeuronal']
                        cluster_str = 'cluster_15CTs'
                    else:
                        CT_list = np.unique(D.ca['cluster_name_'+str(n_cl)+'CTs'])
                        cluster_str = 'cluster_name_'+str(n_cl)+'CTs'
                    CT_list = [ct for ct in CT_list if ct != "none (removed)"]
                    #DF_number = pd.DataFrame(data={"celltype":CT_list.tolist()})
                    DF_number['Number of cells'] = 0
                    #DF_number.set_index("celltype",inplace=True)
                    
                    for ct in CT_list:
                        if ct != "none (removed)":
                            if ct not in DF_number.index.tolist():
                                new_row = pd.Series(data={'up-regulated':0,'down-regulated':0,'total':0,'Number of cells':0},name=ct)
                                #DF_number = DF_number.append(new_row,ignore_index=False) #--> depreciated in version 2.0
                                DF_number = pd.concat([DF_number,pd.DataFrame([new_row])],ignore_index=True)
                            if n_cl==3:
                                n_cells = np.count_nonzero(fnmatch.filter(D.ca[cluster_str],ct.split(" ")[0]+'*'))
                            else:
                                n_cells = np.count_nonzero(fnmatch.filter(D.ca[cluster_str],ct))
                            #matche entry in data frame
                            DF_number.loc[ct,'Number of cells'] = n_cells
                    #lws = [2,1,1]
                    #colors={'total':'black',"down-regulated":'darkred',"up-regulated":'teal'}
                    #ax = DF_number[['Number of cells','total','up-regulated','down-regulated']].plot(x= 'Number of cells', y= ['total','up-regulated','down-regulated'], style = ['o','x','+'], ylabel="Number of DEGs")
                DF_number = get_short_CT_names(DF_number.reset_index(),"celltype")
                DF_number.set_index("celltype_short",inplace = True)
                DF_number.to_csv(path_data+filename_DF_number)

            DF_number = DF_number.reset_index()#.merge(colors, on="celltype").set_index("index")
            DF_number = add_CT_colors_to_DF(DF_number,path_filtered_data,True,n_cl)
            color_dict = DF_number.set_index("celltype_short")["color"].to_dict()
            f, ax = plt.subplots(figsize=(7,4))
            ax.ticklabel_format(style='plain', axis='both')
            sns.regplot(x= 'Number of cells', y= 'total',data=DF_number,fit_reg=True,ax=ax,color="grey")#,hue="celltype",palette = color_dict,ax=ax) 
            sns.scatterplot(x= 'Number of cells', y= 'total',data=DF_number,hue="celltype_short",palette = color_dict,ax=ax) 
            _,xlim_max = ax.get_xlim()
            ax.set(ylabel="Number of DEGs",xlim=[0,xlim_max])
            #ax.set_xlabels("Number of cells")
            #ax = DF_number[['Number of cells','total']].plot(x= 'Number of cells', y= 'total', style = 'o',color="black", xlabel="Number of DEGs",ylabel = "Number of cells",figsize=(3,5),legend=False)
            #ax = remove_frame_keep_axes(ax)
            #_,ylim_max = ax.get_ylim()
            # for i, l in enumerate(ax.lines):
            #     plt.setp(l, linewidth=lws[i])
            #for i,ct_i in enumerate(CT_list):
            #    plt.text(DF_number['Number of cells'][i], ylim_max-0.07*ylim_max, CT_list[i], horizontalalignment='center',verticalalignment='center', fontsize=10, rotation=80)
            
            #add labels:
            #cutoff_number_DEGs_for_label = np.sort(DF_number['total'])[-5]
            #DF_number_sel = DF_number[DF_number['total']>=cutoff_number_DEGs_for_label].copy()
            #for ID, val in enumerate(DF_number_sel['Number of cells']):
            #    print(DF_number_sel.iloc[ID].name)
            #    ax.annotate(DF_number_sel.iloc[ID].name,(DF_number_sel['total'].iloc[ID],val))
    
            #ax.ticklabel_format(useOffset=False,style='plain')
            ax.legend(bbox_to_anchor=(1.37, 1.021))
                
        if opt_save:
            path=path_results+'figures/stacked_bars/'
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)    
            if opt_angel_plot:
                add_str = "_angel"        
            plt.savefig(path+opt_DEG_method+'_'+shrinkage_method+'_num_dereg_genes_'+option+'_'+opt_celltype_groups+'_'+add_figure_str+str(alpha_val).replace('.','_')+add_str+'.pdf', bbox_inches='tight')
        plt.close()
        plt.clf() 

def plot_number_of_celltype_specific_and_shared_pathways(DF_all_pivot,p_val_cutoff,path_results_long,opt_save,add_str,opt_proteomics,opt_proteomics_up_and_down_sep,opt_GOs_only):
    if opt_GOs_only:
        add_str_2 = "_GOs_only"
    else:
        add_str_2 = ""
    loadPltSettings(30,1)
    if opt_proteomics == False and 'proteomics' in DF_all_pivot.columns or opt_proteomics == False and 'proteomics_up' in DF_all_pivot.columns or opt_proteomics == False and 'proteomics_down' in DF_all_pivot.columns:
        #remove proteomics:
        if opt_proteomics_up_and_down_sep:
            DF_all_pivot.drop(columns=['proteomics_up','proteomics_down'],inplace=True)
        else:
            DF_all_pivot.drop(columns=['proteomics'],inplace=True)
    DF_bool_all = DF_all_pivot<=p_val_cutoff
    ct_specific_pws = DF_bool_all.sum(axis=1)[DF_bool_all.sum(axis=1)==1]
    shared_pws = DF_bool_all.sum(axis=1)[DF_bool_all.sum(axis=1)>1]
    n_significant_pws = len(ct_specific_pws)+len(shared_pws)
    if n_significant_pws > 0:
        p_ct_specific_pws = 100*len(ct_specific_pws)/n_significant_pws
        p_shared_pws = 100*len(shared_pws)/n_significant_pws
        plt = loadPltSettings(8,1)
        fig, ax = plt.subplots(figsize=(2,4))
        plt.bar(x=[0.1,0.2],height=[p_ct_specific_pws,p_shared_pws],width=0.05, align='center',tick_label=['cell type specific','shared between cell types'], color = 'grey', edgecolor = 'none')
        ax.set(ylabel='Percentage significant pathways')
        ax = remove_frame_keep_axes(ax)
        #non-scientific y-axis:
        #ax.ticklabel_format(style='plain')
        plt.ticklabel_format(axis='y',style='plain')
        if opt_save:
            path=path_results_long+'figures/stats/'
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            fig.savefig(path + 'n_of_CT_spec_vs_shared_pws_'+add_str.replace(':','_')+add_str_2+'.pdf',bbox_inches='tight')
            plt.close(fig)
            plt.clf()   
    else: 
        print("No plot with the number of cell type specific and shared pathways provided as there were ot significant pathways ("+path_results_long+")!")

#important for heatmap GO term categories on y-axis:
def label_group_bar_table(ax, df, fs, opt_add_line=True,opt_swap_indexes=False,max_label_length = 60, add_offset=0.0, y_axis_reversed=False):
    if opt_swap_indexes:
        idx_range = reversed(range(df.index.nlevels))
    else:
        idx_range = range(df.index.nlevels)
    
        
    xpos = -.15 
    scale = 1./(df.index.size)
    bool_first_level = True
    for level in idx_range:
        pos = df.index.size
        L = label_len(df.index,level)
        if y_axis_reversed==True:
            L.reverse()
        for label, rpos in L:
            #if labl is too long: shorten:
            if len(label)>max_label_length:
                label = label[0:max_label_length]+'...'
                
            pos -= rpos
            
            if bool_first_level:
                # position for the text
                lypos = (pos + .35 * rpos)*scale # pos = 16 and then reduced by 1 in each loop, rpos = 1, rpos = 1/16 
                ax.text(xpos+.1, lypos, label, ha='right', transform=ax.transAxes, fontsize=fs) 
                first_level = level
            else:
                #pathway
                # position for the text
                if opt_add_line:
                    add_line(ax, pos*scale, xpos+.1-add_offset) #xpos-.05
                lypos = (pos + .45 * rpos)*scale # pos = 16 and then reduced by 1 in each loop, rpos = 1, rpos = 1/16 
                ax.text(xpos+.1-add_offset, lypos, label, ha='left', transform=ax.transAxes, fontsize= fs) 
                #print('text position: '+ str(xpos+.1-add_offset))
        bool_first_level=False
        if (level==first_level):
            if (opt_add_line):
                add_line(ax, pos*scale, xpos+.1-add_offset) #xpos-.05
                #print(ax)
                print('first level last line:'+str(xpos+.1-add_offset))

#important for heatmap GO term categories on y-axis:
def add_line(ax, xpos, ypos):
    line = plt.Line2D([ypos, ypos+ .2], [xpos, xpos], color='black', transform=ax.transAxes)
    line.set_clip_on(False)
    ax.add_line(line)

#important for heatmap GO term categories on y-axis:
def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def plot_heatmap_pathways_from_DEGs_per_celltype(DF_sign,n_clusters, opt_save, folder, variables_groups, n_genes_names, data_path, opt_transcriptomics_up_and_down_sep):
    #for up and down seperately:
    if opt_transcriptomics_up_and_down_sep:
        for reg_str in ['up','down']:
            if reg_str=='up':
                colors = ['#FFECEC','#8B0000']
            else:
                colors = ['#DCF6F6','#008080']
                
            #DF_sign.sort_values(by=variables_groups[1],axis=0,inplace=True)
            
            #drop the rows not interested in --> ending with '_'+reg_str+'_sig'
            DF_sign = DF_sign[DF_sign['Celltype'].str.endswith('_'+reg_str+'_sig')]
            if len(DF_sign)==0:
                continue
            else:
                DF_sign[n_genes_names]=DF_sign[n_genes_names].astype('int')
                DF_sign['genes.in.geneset']=DF_sign['genes.in.geneset'].astype('int')
                DF_sign['percentage_genes_overlap'] = np.round(DF_sign[n_genes_names]/DF_sign['genes.in.geneset']*100,0).astype('int')
            
                for normalization in ['','geneset_size']:
                    for aggr in ['sum','mean']:
                        if normalization=='geneset_size' and aggr == 'mean':
                            heatmap_data = pd.pivot_table(DF_sign, values = 'percentage_genes_overlap', index=variables_groups, columns='Celltype', aggfunc = aggr)
                            add_str = 'percentage_'
                        else:
                            heatmap_data = pd.pivot_table(DF_sign, values = n_genes_names, index=variables_groups, columns='Celltype', aggfunc = aggr)
                            add_str = 'number_'

                        #remove lines which are exclusively nan
                        heatmap_data=heatmap_data.dropna(how='all')
                        #resort columns in heatmap_data:
                        CT_order = get_CT_order(n_clusters,data_path,True)
                        if "Patient_clustering/" in folder:
                            CT_order =  [string1+'__'+reg_str+'_sig' for string1 in CT_order]
                        else:
                            CT_order =  [string1.replace('-','_').replace(' ','_')+'_'+reg_str+'_sig' for string1 in CT_order]
                        
                        #add missiong columns:
                        if len(CT_order)>np.shape(heatmap_data)[1]:
                            #add missing CT as column with nans
                            columns_missing = [ct for ct in CT_order if ct not in heatmap_data.columns.tolist()]   
                            for col_m in columns_missing:
                                heatmap_data[col_m]=np.nan
                        #sort columns according to CT order:
                        columns_ordered = [ct for ct in CT_order if ct in heatmap_data.columns.tolist()]    
                        heatmap_data = heatmap_data[columns_ordered]
                        
                        #sort values by group
                        heatmap_data.sort_values(by=variables_groups[1],axis=0,inplace=True)
                        if np.shape(heatmap_data)[0]>0:
                            #corresponding list of GO term info
                            plt = loadPltSettings(25,1)
                            fig = plt.figure(figsize=(80,80))
                            ax = fig.add_subplot(111)
                            formatter = tkr.ScalarFormatter(useMathText=True)
                            formatter.set_scientific(False)
                            #later: formatter.set_powerlimits((0, overall max))
                            #more fine tuning for legend: https://stackoverflow.com/questions/43697234/making-colorbar-values-integer-in-a-heatmap-matplotlib-seaborn
                            #if subfolder == 'alpha30':
                            #define color bar label str:
                            if aggr=='sum':
                                cb_label_str = "Total number "
                            elif aggr == 'mean' and normalization=='geneset_size':
                                cb_label_str = " Average percentage "
                            else:
                                cb_label_str = "Average number "
                            cb_label_str = cb_label_str + "of " + reg_str+"-regulated genes"
                            if min(np.shape(heatmap_data))>0:
                                sns.heatmap(heatmap_data, 
                                            vmin=0, 
                                            annot=True, 
                                            cmap=get_continuous_cmap(colors), 
                                            linewidths=.5, 
                                            linecolor='black', 
                                            clip_on=False,
                                            fmt='g',
                                            cbar_kws={"format": formatter, "label": cb_label_str}, 
                                            square=True,
                                            annot_kws={"size":12})
    
                            #Below 3 lines remove default labels
                            labels = ['' for item in ax.get_yticklabels()]
                            ax.set_yticklabels(labels, fontsize=12)
                            ax.set_ylabel('')
                            x_labels = [item.get_text().replace('__'+reg_str+'_sig','') for item in ax.get_xticklabels()]
                            ax.set_xticklabels(x_labels, fontsize=25)
                            label_group_bar_table(ax, heatmap_data,fs=30,add_offset=0.4)
                            fig.subplots_adjust(bottom=.1*heatmap_data.index.nlevels)
                            if opt_save:
                                path=folder+'figures/heatmaps/'
                                isExist = os.path.exists(path)
                                if not isExist:
                                    os.makedirs(path)
                                os.chdir(path)
                                fig.savefig(path+'heatmap_'+add_str+reg_str+'_regulated_'+aggr+'.pdf',bbox_inches='tight')
                                plt.close(fig)
                                plt.clf() 
    else:
        colors = ['#E4E4E4','#040404']
        if len(DF_sign)>0:
            DF_sign[n_genes_names]=DF_sign[n_genes_names].astype('int')
            DF_sign['genes.in.geneset']=DF_sign['genes.in.geneset'].astype('int')
            DF_sign['percentage_genes_overlap'] = np.round(DF_sign[n_genes_names]/DF_sign['genes.in.geneset']*100,0).astype('int')
        
            for normalization in ['','geneset_size']:
                for aggr in ['sum','mean']:
                    if normalization=='geneset_size' and aggr == 'mean':
                        heatmap_data = pd.pivot_table(DF_sign, values = 'percentage_genes_overlap', index=variables_groups, columns='Celltype', aggfunc = aggr)
                        add_str = 'percentage_'
                    else:
                        heatmap_data = pd.pivot_table(DF_sign, values = n_genes_names, index=variables_groups, columns='Celltype', aggfunc = aggr)
                        add_str = 'number_'
                        
                        
                    # TO DO:
                        #cell type names (columns) in heatmap_data have different ending as cell types in CT_order
                        
                    #remove lines which are exclusively nan
                    heatmap_data=heatmap_data.dropna(how='all')
                    #resort columns in heatmap_data:
                    CT_order = get_CT_order(n_clusters,data_path,True)
                    CT_order =  [string1.replace('-','_').replace(' ','_')+'_sig' for string1 in CT_order]
                    
                    #add missiong columns:
                    if len(CT_order)>np.shape(heatmap_data)[1]:
                        #add missing CT as column with nans
                        columns_missing = [ct for ct in CT_order if ct not in heatmap_data.columns.tolist()]   
                        for col_m in columns_missing:
                            heatmap_data[col_m]=np.nan
                    #sort columns according to CT order:
                    columns_ordered = [ct for ct in CT_order if ct in heatmap_data.columns.tolist()]    
                    heatmap_data = heatmap_data[columns_ordered]
                    
                    #sort values by group
                    heatmap_data.sort_values(by=variables_groups[1],axis=0,inplace=True)
                    if np.shape(heatmap_data)[0]>0:
                        #corresponding list of GO term info
                        plt = loadPltSettings(50,1)
                        fig = plt.figure(figsize=(80,60))
                        ax = fig.add_subplot(111)
                        formatter = tkr.ScalarFormatter(useMathText=True)
                        formatter.set_scientific(False)
                        #later: formatter.set_powerlimits((0, overall max))
                        #more fine tuning for legend: https://stackoverflow.com/questions/43697234/making-colorbar-values-integer-in-a-heatmap-matplotlib-seaborn
                        #if subfolder == 'alpha30':
                        #define color bar label str:
                        if aggr=='sum':
                            cb_label_str = "Total number "
                        elif aggr == 'mean' and normalization=='geneset_size':
                            cb_label_str = " Average percentage "
                        else:
                            cb_label_str = "Average number "
                        cb_label_str = cb_label_str + "of deregulated genes"
                        if min(np.shape(heatmap_data))>0:
                            sns.heatmap(heatmap_data, vmin=0, annot=True, cmap=get_continuous_cmap(colors), linewidths=.5, linecolor='black', clip_on=False,fmt='g',cbar_kws={"format": formatter, "label": cb_label_str}, square=True)

                        #Below 3 lines remove default labels
                        labels = ['' for item in ax.get_yticklabels()]
                        ax.set_yticklabels(labels)
                        ax.set_ylabel('')
                        label_group_bar_table(ax, heatmap_data,fs=40,add_offset=0.2)
                        fig.subplots_adjust(bottom=.1*heatmap_data.index.nlevels)
                        if opt_save:
                            path=folder+'figures/heatmaps/'
                            isExist = os.path.exists(path)
                            if not isExist:
                                os.makedirs(path)
                            os.chdir(path)
                            fig.savefig(path+'heatmap_'+add_str+'_deregulated_'+aggr+'.pdf',bbox_inches='tight')
                            plt.close(fig)
                            plt.clf() 

def plot_PCA_of_pseudobulk_data_for_gene_selection(cts,options, path_info, n_cl, factor, genelist_kind, add_fig_str):

    #initialize figures
    plt = loadPltSettings(FS,10)
    plt.rcParams['axes.labelsize']= 20
    plt.ticklabel_format(style='plain') 
    n_cols = int(np.ceil(n_cl/4))

    f_expVar, ax_expVar = plt.subplots(4, n_cols, figsize=(35,30))
    f_expVar.subplots_adjust(hspace = .35, wspace = 0.15) 

    f_pc1_2, ax_pc1_2 = plt.subplots(4, n_cols, figsize=(35,30))
    f_pc1_2.subplots_adjust(hspace = .35, wspace = 0.15) 
    
    r_id = 0
    c_id = 0
    
    #get list of samples and their group:
    filename_sample_info = path_info["path_project"]+'data/Sample_information/sample_info_SUN_SCZ.xlsx'
    DF_sample_info = pd.read_excel(filename_sample_info,engine="openpyxl")
    DF_group_sample = DF_sample_info[['donor_ID_internal_8','scz_status']]
    DF_group_sample.rename(columns={'donor_ID_internal_8':'Sample_ID','scz_status':'Group'},inplace=True)
    DF_group_sample.dropna(axis = 'index',inplace=True)
    #rename columns to Group and Sample_ID
        
    for ct_id,ct in enumerate(cts):
        #get pseudobulk data for current ct:
        filename_aggr_counts_ct_i = path_info["path_project"]+'/data/CT_clustered_loom_formatted_data_cellranger/'+str(n_cl)+'_CTs/Aggregated_data/Aggregated_counts_cellranger_pagoda_TH_and_D_adj_filtered_conos_cluster_based_'+ct+'_sum.csv'
        aggr_counts_ct_i = pd.read_csv(filename_aggr_counts_ct_i)
        if add_fig_str == "_mean":
            #normalize aggregated data for each donor by its total sum of counts --> leads to mean instead of sum aggregation
            sample_columns = aggr_counts_ct_i.columns[aggr_counts_ct_i.columns.str.startswith('S')]
            aggr_counts_ct_i[sample_columns] = aggr_counts_ct_i[sample_columns]/aggr_counts_ct_i[sample_columns].sum()
        #create lists of genes for current factor
        if genelist_kind=="input_GSA_of_factors":
            genes = []
            if options["opt_up_and_down_sep"]==True:
                reg_strings = ['up','down'] 
            else:
                reg_strings = ['up_or_down'] 
            for reg_str in reg_strings:
                file_gene_list = path_info["path_results_subdirectory"] + 'gene_list_' + ct + '_' + reg_str + '.csv'
                genes_tmp = pd.read_csv(file_gene_list,header=None,names=['Gene'])['Gene'].tolist()
                genes = genes + genes_tmp
        elif genelist_kind == "genes_in_loadings_matrix_of_factor":
            if ct_id==0: # same for all cell types, read only once
                file_gene_list = path_info["path_results_subdirectory"] + 'factor_' + factor + '_loadings_matrix.csv'
                genes_tmp = pd.read_csv(file_gene_list,header="infer")
                genes = genes_tmp["Unnamed: 0"].tolist()
        ax_expVar,ax_pc1_2,r_id,c_id = plot_PCA_of_pseudocounts_for_ct_and_selected_genes(aggr_counts_ct_i,ct_id,ct,cts,genes,DF_group_sample,ax_expVar,ax_pc1_2,r_id,c_id,n_cols)
        if ct_id+1==len(cts):
            #last cell type
            if r_id<4:
                #delete axes of remaining subplots
                for c_id_del in range(c_id,n_cols-1,1):
                    ax_pc1_2[r_id,c_id_del].set_axis_off()
                    ax_expVar[r_id,c_id_del].set_axis_off()

            #create path
            if genelist_kind=="input_GSA_of_factors":
                path = path_info["path_results_subdirectory"]+"figures/PCA/"
            elif genelist_kind == "genes_in_loadings_matrix_of_factor":
                path = path_info["path_results_subdirectory"]+"PCA/"
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)

            #save figures:
            f_expVar.suptitle('Explained variance of PCA of pseudobulk counts filtered for '+factor+'-genes')
            f_expVar.savefig(path+'expl_var_all_cts_factor_'+factor+"_genes"+add_fig_str+'_'+genelist_kind+'.pdf', bbox_inches='tight')
            plt.close(f_expVar)

            f_pc1_2.suptitle('PCA of pseudobulk counts filtered for '+factor+'-genes')
            f_pc1_2.savefig(path+'PC1_vs_PC2_all_cts_factor_'+factor+"_genes"+add_fig_str+'_'+genelist_kind+'.pdf', bbox_inches='tight')
            #plt.savefig(path_results_figures+'PCA_'+ca_ann[k]+'.png', bbox_inches='tight')
            plt.close(f_pc1_2)

def scatter_text(x_var, y_var, text_column, hue_variable, dataset, title, xlabel, ylabel, axes, opt_label_outliers_only,dot_size=12, text_size='large',font_color="black"):
    #Scatter plot with country codes on the x y coordinates
    #  Based on this answer: https://stackoverflow.com/a/54789170/2641825
    # Create the scatter plot
    p1 = sns.scatterplot(data=dataset, x=x_var, y=y_var, hue=hue_variable, s = dot_size, palette=['steelblue','darkorange'],ax=axes)
    # Add text besides each point
    #add_option only for outliers:
    for line in range(0,dataset.shape[0]):
         if opt_label_outliers_only==True:
            if (dataset[x_var][line]-dataset[x_var].mean())**2/(4*dataset[x_var].std()) + (dataset[y_var][line]-dataset[y_var].mean())**2/(4*dataset[y_var].std()) > 1:
                p1.text(dataset[x_var][line]+0.01*(np.max(dataset[x_var])-np.min(dataset[x_var])), dataset[y_var][line], 
                        dataset[text_column][line], horizontalalignment='left', 
                        size=text_size, color=font_color)#, weight='semibold')
         else:
            p1.text(dataset[x_var][line]+0.01*(np.max(dataset[x_var])-np.min(dataset[x_var])), dataset[y_var][line], 
                    dataset[text_column][line], horizontalalignment='left', 
                    size=text_size, color=font_color)#, weight='semibold')
    # Set title and axis labels
    axes.set_title(title,fontsize=20)
    axes.set_xlabel(xlabel,fontsize=20)
    axes.set_ylabel(ylabel,fontsize=20)
    return p1

def plot_PCA_of_pseudocounts_for_ct_and_selected_genes(aggr_counts_ct_i,ct_id,ct,cts,genes,DF_group_sample,ax_expVar,ax_pc1_2,r_id,c_id,n_cols):
        #subset pseudobulk data
        aggr_counts_ct_i_sel = aggr_counts_ct_i[aggr_counts_ct_i['Accession'].isin(genes)]
        aggr_counts_ct_i_sel.drop(labels=["Accession","Gene"],axis='columns',inplace=True)
        
        #plot pca
        pca = PCA(n_components=20)
        pca.fit(np.transpose(aggr_counts_ct_i_sel))
        print("Variance explained by all 20 principal components = ", sum(pca.explained_variance_ratio_*100))
        
        ax_expVar[r_id, c_id].plot(np.cumsum(pca.explained_variance_ratio_))
        ax_expVar[r_id, c_id].set_title(ct)
        ax_expVar[r_id, c_id].ticklabel_format(useOffset=False, style='plain')

        aggr_counts_ct_i_sel_pca=pca.transform(np.transpose(aggr_counts_ct_i_sel))

        PCA_R_tmp = pd.DataFrame(data = {'Sample_ID': aggr_counts_ct_i_sel.columns.tolist(), 'PC1': aggr_counts_ct_i_sel_pca[:,0], 'PC2': aggr_counts_ct_i_sel_pca[:,1]})
        #To Do:  add Group to dataframe
        PCA_R = pd.merge(left = PCA_R_tmp,right=DF_group_sample, how="inner", on="Sample_ID")

        if r_id==3:
            ax_expVar[r_id, c_id].set_xlabel('Number of components')
        else:
            ax_expVar[r_id, c_id].set_xlabel('')

        if c_id==0: 
            ax_expVar[r_id, c_id].set_ylabel('Explained variance')
        else:
            ax_expVar[r_id, c_id].set_ylabel('')


        scatter_text(x_var='PC1',
                     y_var='PC2',
                     text_column='Sample_ID',
                     hue_variable='Group',
                     dataset=PCA_R, 
                     title=ct,
                     xlabel='PC 1 ('+str(np.round(pca.explained_variance_ratio_[0]*1000)/10)+' %)',
                     ylabel='PC 2 ('+str(np.round(pca.explained_variance_ratio_[1]*1000)/10)+' %)',
                     axes=ax_pc1_2[r_id,c_id],
                     opt_label_outliers_only = False,
                     dot_size=20,
                     text_size='small',
                     font_color = 'lightgrey')

        if ct_id+1==len(cts):
            #last cell type:
            ax_pc1_2[r_id,c_id].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,scatterpoints=1, fontsize=20)
        else:
            ax_pc1_2[r_id,c_id].legend_ = None

        #update column and row ids
        c_id = c_id+1
        if np.mod(c_id,n_cols)==0:
            r_id = r_id+1
            c_id = 0

        return ax_expVar,ax_pc1_2,r_id,c_id

def plot_detailed_stacked_bar_DEGs(DF_DEG_results_detailed,opt_save, path_results,opt_DEG_method,shrinkage_method, option,opt_celltype_groups,add_figure_str, alpha_val,n_clusters,path_filtered_data):
    DF_DEG_results_detailed = DF_DEG_results_detailed.reset_index(drop=False).sort_values(by="total",ascending=False)
    if "celltype_short" in DF_DEG_results_detailed.columns:
        DF_DEG_results_detailed.set_index(keys="celltype_short",inplace=True)
    elif "celltype_names_short" in DF_DEG_results_detailed.columns:
        DF_DEG_results_detailed.set_index(keys="celltype_names_short",inplace=True)
    DF_DEG_results_detailed.index = DF_DEG_results_detailed.index.str.replace('_',' ')
    filename_fig = opt_DEG_method+'_'+shrinkage_method+'_num_dereg_genes_'+option+'_'+opt_celltype_groups+'_'+add_figure_str+str(alpha_val).replace('.','_')+'_gene_type.pdf'
    plot_detailed_stacked_bar(DF_DEG_results_detailed,opt_save,path_results,filename_fig)

def plot_detailed_stacked_bar(DF_results_detailed,opt_save,path_results,filename_fig):
    axes=DF_results_detailed[['down-regulated protein coding','down-regulated lncRNA','up-regulated protein coding','up-regulated lncRNA']].plot(kind='barh',
                                                                                                                                                stacked=True,
                                                                                                                                                xlabel='Number of genes',
                                                                                                                                                color={'down-regulated protein coding':'teal','down-regulated lncRNA':'#104E4C','up-regulated protein coding':'darkred','up-regulated lncRNA':'#5B0B0B'},
                                                                                                                                                edgecolor="white")
    axes.ticklabel_format(style='plain', axis='x') 
    axes = remove_frame_keep_axes(axes)
    plt.gca().invert_yaxis()

    if opt_save:
        path=path_results+'figures/stacked_bars/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)            
        plt.savefig(path+filename_fig, bbox_inches='tight')
        plt.close()
        plt.clf() 

def plot_DEGs_CT_overlap_as_upset(DF_DESeq2,opt_DEG_method,opt_save,path_results,aggregation_method,add_str,alpha_val):
    #get short CT names
    DF_DESeq2 = get_short_CT_names(DF_DESeq2,"celltype_names")
    if alpha_val < 0.2:
        min_subset_size=3
    else:
        min_subset_size=7
    if np.shape(DF_DESeq2)[0]>0:
        if opt_DEG_method == 'Wilcoxon':
            columns = ['Cell type', 'names']
            #add deregulation column
            DF_DESeq2.loc[(DF_DESeq2.padj>alpha_val),'deregulation']='none'
            DF_DESeq2.loc[(DF_DESeq2.padj<=alpha_val) & (DF_DESeq2.logfoldchanges>0),'deregulation']='up'
            DF_DESeq2.loc[(DF_DESeq2.padj<=alpha_val) & (DF_DESeq2.logfoldchanges<0),'deregulation']='down'
        elif opt_DEG_method=='DESeq2':
            if 'shrinkage' in DF_DESeq2.columns:
                if 'celltype_names_short' in DF_DESeq2.columns:
                    columns = ['celltype_names_short', 'genes']
                elif 'celltype_name_short' in DF_DESeq2.columns:
                    columns = ['celltype_name_short', 'genes']
                elif 'celltype_name' in DF_DESeq2.columns:
                    columns = ['celltype_name', 'genes']
                elif 'celltype_names' in DF_DESeq2.columns:
                    columns = ['celltype_names', 'genes']
                else:
                    columns = ['celltype', 'genes']
            else:
                if 'celltype_names_short' in DF_DESeq2.columns:
                    columns = ['celltype_names_short', 'Gene']
                else:
                    columns = ['celltype_names', 'Gene']
            #add deregulation column
            DF_DESeq2.loc[(DF_DESeq2.padj>alpha_val),'deregulation']='none'
            DF_DESeq2.loc[(DF_DESeq2.padj<=alpha_val) & (DF_DESeq2.log2FoldChange>0),'deregulation']='up'
            DF_DESeq2.loc[(DF_DESeq2.padj<=alpha_val) & (DF_DESeq2.log2FoldChange<0),'deregulation']='down'
            
            DF_DESeq2[columns[0]] = DF_DESeq2[columns[0]].str.replace("_"," ")
        else: 
            print("DEG_method "+opt_DEG_method+" is not specified!")
        
        reg_str = ['up','down','all_seperately','all_together']
        colors = ['darkred','teal','black','black']
        for i,r_str in enumerate(reg_str):
            if r_str=='all_seperately':
                DF_DESeq2[columns[0]] = DF_DESeq2[columns[0]]+ ' ' + DF_DESeq2['deregulation']
                DF_sel = DF_DESeq2.copy()
            elif r_str=='all_together':
                DF_sel = DF_DESeq2.copy()
                #put up and down together:
                CTs = DF_sel['celltype_names'].unique().tolist()
                CTs_raw = np.unique([item.split(' up', 1)[0] for item in CTs if item.endswith(' up')] + [item.split(' down', 1)[0] for item in CTs if item.endswith(' down')]).tolist()
                #sum up up- and down-reg:
                for ct in CTs_raw:
                    if ct+' up' in CTs:
                        DF_sel.loc[DF_sel.celltype_names==ct+' up','celltype_names'] = ct
                    if ct+' down' in CTs:
                        DF_sel.loc[DF_sel.celltype_names==ct+' down','celltype_names'] = ct
                CTs=CTs_raw     
            else:
                DF_sel = DF_DESeq2[DF_DESeq2['deregulation']==r_str][columns].copy()
            if np.shape(DF_sel)[0]>0:
                CTs = np.unique(DF_sel[columns[0]].tolist()).tolist()
                for ct_i in CTs:
                    DF_sel[ct_i]=False
                    DF_sel.loc[DF_sel[columns[0]]==ct_i,ct_i]=True
                DF_sel.drop(columns = columns[0],inplace=True)
                DF_sel.reset_index(inplace=True)
                DF_sel.drop(columns = 'index',inplace=True)
                #sum rows if same gene:
                DF_sel_red = DF_sel.groupby(columns[1])[CTs].agg('sum')
                #transform to bool:
                DF_sel_red = DF_sel_red==1
                #rename columns (add spaces)
                col_names=DF_sel_red.columns.tolist()
                #for col in col_names:
                #    eval("DF_sel_red.rename(columns={\'"+col+"\':\'.               "+col+"\'},inplace=True)")
                #DF_sel_red.rename(columns={'Excitatory Layer 2-3 LINC00507 RPL9P17 neurons II':'.               Excitatory Layer 2-3 LINC00507 RPL9P17 neurons II'})
                #plot:
                params = {'font.size':10, 'axes.labelsize':10, 'xtick.labelsize':10, 'ytick.labelsize':10}
                CTs_updated = DF_sel_red.columns.tolist()
                with plt.rc_context(params):
                    plot(from_indicators(CTs_updated,data=DF_sel_red), show_counts=True, sort_by = 'cardinality', sort_categories_by='cardinality', min_subset_size=min_subset_size, facecolor=colors[i])
                    
                if opt_save:
                    path=path_results+'figures/upset/'
                    isExist = os.path.exists(path)
                    if not isExist:
                        os.makedirs(path)
                    plt.savefig(path+'upset_'+r_str+'_regulated_DEGs_CT_overlap'+opt_DEG_method+'_'+aggregation_method+add_str+'_'+str(alpha_val).replace('.','_')+'.pdf', bbox_inches='tight')
                plt.close()
                plt.clf() 
        print(DF_sel)
        #To do: plot number of DEGs for DEGs shared across cell types up-reg, shared across cell types down-reg, shared across cell types up- and down-reg, and cell type specific
        #Group by two columns and count the occurrences of each combination:

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

def plot_number_pathways_per_gene(DF_sign_pw_per_gene, opt_save, folder,n_clusters,path_colors,add_str):
    if np.shape(DF_sign_pw_per_gene)[0]>0:
        plt = loadPltSettings(12,1)
        plt.rcParams['axes.formatter.useoffset'] = False
        fig, ax = plt.subplots(figsize=(9,6))
        #plt.gca().set_adjustable("box")
        DF_sign_pw_per_gene['celltype'] = DF_sign_pw_per_gene['celltype'].str.replace('_',' ').str.replace('2 3','2-3').str.replace('3 4','3-4').str.replace('3 6','3-6').str.replace('5 6','5-6')
        DF_sign_pws_per_gene_up = DF_sign_pw_per_gene[DF_sign_pw_per_gene['regulation']=='up'].pivot_table(values='Number_significant_pathways_per_gene',index='celltype',aggfunc='mean')
        DF_sign_pws_per_gene_down = DF_sign_pw_per_gene[DF_sign_pw_per_gene['regulation']=='down'].pivot_table(values='Number_significant_pathways_per_gene',index='celltype',aggfunc='mean')
        DF_sign_pws_per_gene_up = DF_sign_pws_per_gene_up.rename(columns={'Number_significant_pathways_per_gene': 'Number_significant_pathways_per_up-regulated_gene'})
        DF_sign_pws_per_gene_down = DF_sign_pws_per_gene_down.rename(columns={'Number_significant_pathways_per_gene': 'Number_significant_pathways_per_down-regulated_gene'})
        DF = pd.concat([DF_sign_pws_per_gene_up,DF_sign_pws_per_gene_down], axis=1)
        CT_color_df = get_CT_colors_as_df(n_clusters,path_colors)
        CT_color_df = CT_color_df.set_index('celltype')
        df = pd.concat([DF, CT_color_df], axis=1)
        df.reset_index(inplace=True)
        df = get_short_CT_names(df,"celltype")
        df.set_index("celltype_short",inplace=True)
        #df = pd.merge(DF, CT_color_df, how='left', left_on='celltype', right_index=True)
        if df['Number_significant_pathways_per_up-regulated_gene'].sum()>0 or df['Number_significant_pathways_per_down-regulated_gene'].sum()>0:
            for ticker,row in df.iterrows():
                ax.scatter(x=row['Number_significant_pathways_per_up-regulated_gene'],y=row['Number_significant_pathways_per_down-regulated_gene'],color=row['color'],label = ticker)
            add_identity(ax, color='black', ls='--')
            #ax.axis('equal')
            _,y2 = ax.get_ylim()
            _,x2 = ax.get_xlim()
            ax.set_xlim([0,int(math.ceil(np.max([x2,y2])))])
            ax.set_ylim([0,int(math.ceil(np.max([x2,y2])))])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.set_xlabel('Number significant pathways per up-regulated gene')
            ax.set_ylabel('Number significant pathways per down-regulated gene')
            ax.set_aspect('equal', 'box') #ax.axis('equal')
            _,b=ax.get_xlim()
            _,d=ax.get_ylim()
            m = np.max([b,d])
            ax.set_ylim([0,m])
            ax.set_xlim([0,m])
            ax.ticklabel_format(axis='both',style='plain')
            ax = remove_frame_keep_axes(ax)
            if opt_save:
                path=folder+'figures/stats/'
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                fig.savefig(path + 'number_sign_pathways_per_gene_'+add_str.replace(':','_')+'.pdf',bbox_inches='tight')
                plt.close(fig)
                plt.clf()       

def plot_comparison_with_randomized_group_labels(DF_true,path_results,opt_save):
    DF_true.sort_values('log2FoldChange_abs',axis=0,ascending=False,inplace=True)
    DF_true['Rank'] = range(1,np.shape(DF_true)[0]+1)
    plt = loadPltSettings(14,3)
    ax = DF_true.plot(x = 'Rank',y='P_random_DEG', kind='scatter')
    ax.set_xlabel('Rank of abs. log2(FC) from high to low of DEGs from original data')
    ax.set_ylabel('% of permuted datasets where a gene is wrongly identified as DEG')
    if opt_save:
        plt.savefig(path_results+'Comparison_wíth_randomized_group_labels.png',dpi=150,bbox_inches='tight')
        plt.close()
        plt.clf() 

def plot_percentage_wrongly_assigned_group_labels(path_results_short,opt_save): 
    P = get_percentage_wrongly_assigned_group_labels(path_results_short+'randomized_group_labels/')
    plt = loadPltSettings(14,3)
    plt.hist(P)
    plt.xlabel("Percentage of wrongly assigned group labels")
    plt.ylabel("Frequency")
    if opt_save:
        plt.savefig(path_results_short+'randomized_group_labels/'+'Percentage_wrongly_assigned_group_labels.png',dpi=150,bbox_inches='tight')
        plt.close()
        plt.clf()

def plot_venn_geneMatrix(DF,genelist,genelist_name,opt_save,path_results,opt_analysis_method,aggregation_method,shrinkage_method,wilcoxon_correction_method,alpha_val):
    if opt_analysis_method=='DESeq2':
        DF = DF[DF['shrinkage']==shrinkage_method].copy()
        s1=DF['genes']
    else:
        s1=DF['names']
    plt.figure(figsize=(4,4))
    out = venn2([set(s1),set(genelist)],
          set_colors=('#3E64AF', '#3EAF5D'),
          set_labels=[opt_analysis_method+' '+ wilcoxon_correction_method,genelist_name])
    venn2_circles([set(s1), set(genelist)], lw=0.7)
    for text in out.set_labels:
        text.set_fontsize(10)
    for x in range(len(out.subset_labels)):
        if out.subset_labels[x] is not None:
            out.subset_labels[x].set_fontsize(10)
    plt.title('Overlab of '+opt_analysis_method+' '+ wilcoxon_correction_method +' with '+genelist_name)
    #plt.show()  
    if opt_save:
        path=path_results+'venn_genelists/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'venn_'+opt_analysis_method+'_'+wilcoxon_correction_method+aggregation_method+'_'+shrinkage_method+'_with_'+genelist_name+'_'+str(alpha_val)+'.pdf', bbox_inches='tight')
    plt.close()
    plt.clf()  
      
def plot_correlation_DEG_pvalues_with_genelist_p_or_qvalues(DF_DEG,genelist,genelist_name,p_or_qvalues,opt_save,path_results,CT_i,alpha_val,cutoff_q_p):
    if genelist_name in ['asd.wes','scz.wes']:
        string1 = 'q-value '
    else:
        string1 = 'p-value '
    min_number_pairs_to_plot_corr = 20
    #check if qvalues empty:
    if len(p_or_qvalues)>0:
        #filter DF for genes in genelist
        #make sure order of genes is the same
        #best way: create another dataframe and merge them 
        DF = get_genelist_DEG_DF(DF_DEG, genelist, genelist_name, p_or_qvalues, string1)
        #change color for p_val < 0.5 and > 0.05 and q-val < 0.5 to black:
        DF.loc[(DF['padj']>alpha_val) & (DF['padj']<cutoff_q_p),'Color_sign'] = 'black'
        DF.loc[(DF[string1+genelist_name]<cutoff_q_p) & (DF['padj']>alpha_val), 'Color_sign'] = 'black'
        # DF.loc[(DF['padj']>alpha_val) & (DF['padj']<cutoff_q_p) & (DF[string1+genelist_name]<cutoff_q_p),'Color_sign'] = 'black'
        #treat up- and downregulated genes separately
        reg_str = ['up','down']
        for i,r in enumerate(reg_str):
            if r=='up':
                DF_r = DF[DF['log2FoldChange']>0].copy()
            else:
                DF_r = DF[DF['log2FoldChange']<0].copy()
            
            #plot qvalues vs  p values
            if len(p_or_qvalues)>min_number_pairs_to_plot_corr:
                ax = DF_r.plot(x='-log10(padj)', y='-log10('+string1+genelist_name+')', kind= 'scatter', c = 'Color_sign', edgecolors=[0, 0, 0, 0], s=3)
                ax.set_xlabel('-log10(adjusted p-values DEG analysis)')
                ax.set_xlim([0,1])
                ax.set_ylim([-0.15,1.05])
                ax.invert_xaxis()
                ax.invert_yaxis()
                #label genes with sign pvalues or qvalues with gene name
                for p, q in DF_r[(DF_r['Color_sign']!='grey') & (DF_r['Color_sign']!='black')][['padj',string1+genelist_name]].iterrows():
                    ax.annotate(p, q, fontsize=8)
            #calculate correlation
            # DF_sel = DF_r[DF_r['Color_sign']!='grey'].copy()
            if r=='up':
                n_datapoints_up = len(DF_r)
                corr_spearman_up = np.corrcoef(x=DF_r['-log10(padj)'],y=DF_r['-log10('+string1+genelist_name+')'])[0,1]
                if len(p_or_qvalues)>min_number_pairs_to_plot_corr:
                    ax.text(0.75,-0.1,"rho = "+str(np.round(corr_spearman_up,3)),fontsize=12)
            else:
                n_datapoints_down = len(DF_r)
                corr_spearman_down = np.corrcoef(x=DF_r['-log10(padj)'],y=DF_r['-log10('+string1+genelist_name+')'])[0,1]
                if len(p_or_qvalues)>min_number_pairs_to_plot_corr:
                    ax.text(0.75,-0.1,"rho = "+str(np.round(corr_spearman_down,3)),fontsize=12)
                #perform correlation test and add result as text
            if len(p_or_qvalues)>min_number_pairs_to_plot_corr and opt_save:
                path=path_results+'figures/corr_with_lists/'
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                plt.savefig(path+CT_i.replace('-','_')+'_'+r+'_vs_'+genelist_name.replace('.','_')+'.pdf', bbox_inches='tight')
                plt.close()
                plt.clf()  
            
            plot_violin_pvals_genelists(DF_r,r,genelist_name,alpha_val,CT_i,opt_save,path_results,string1)
        
        return(corr_spearman_up,corr_spearman_down,n_datapoints_up, n_datapoints_down)

def plot_violin_pvals_genelists(DF,reg_str,genelist_name,alpha_val,CT_i,opt_save,path_results,string1):
    DF['Group'] = ['False']*np.shape(DF)[0]
    DF.loc[DF[string1+genelist_name]<alpha_val,'Group']=True
    DF_s = DF[DF['-log10(padj)']<=alpha_val].copy()
    if DF_s.empty==False:
        if reg_str == 'up':
            c = 'darkred' 
        else:
            c = 'teal'
        #[x1,x2] = ax1.get_xlim()
        ax1 = sns.swarmplot(x="Group",y="-log10(padj)",data=DF_s,color=c)
        #ax1.set_xlim([x1,x2])
        sns.violinplot(x="Group",y="-log10(padj)",data=DF,title=genelist_name,ax=ax1,cut=0)
    else:
        ax1 = sns.violinplot(x="Group",y="-log10(padj)",data=DF,title=genelist_name,cut=0)
    #plt.setp(ax1.collections, alpha=.3)
    for art in ax1.get_children():
        if isinstance(art, PolyCollection):
            art.set_alpha(0.3)

    #for g, p in DF_s[['Group','padj']].iterrows():
    #    ax1.annotate(g, p, fontsize=8)
    ax1.invert_xaxis()
    ax1.set_ylim([0,1])
    #save plot (optionally)
    if opt_save:
        path=path_results+'figures/violin_genelists/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'violin_'+CT_i+'_'+reg_str+'_DEG_pvals_for_'+genelist_name+'.pdf', bbox_inches='tight')
    plt.close()
    plt.clf()  

def plot_volcano(df_DESeq2_all_genes,CT_i,path_results,alpha_val,opt_save,annotations="very_small_pvals"):
    #plt = loadPltSettings(12,4)
    #fig, ax = plt.subplots()
    if annotations=="very_small_pvals" or annotations=="small_pvals":
        df_DESeq2_all_genes.loc[df_DESeq2_all_genes['padj']>alpha_val,"Color_sign"] = "lightgrey"
        ax = df_DESeq2_all_genes.plot(kind= 'scatter', x='log2FoldChange', y='negative_log10_padj', c="Color_sign",edgecolor='none', s=7, title = CT_i.replace('_',' '))
        #ax = df_DESeq2_all_genes.plot(kind= 'scatter', x='log2FoldChange', y='negative_log10_padj', c="Color_sign",edgecolor='none', s=7, title = CT_i.replace('_',' '))
    else:
        df_DESeq2_all_genes[annotations + '_color'] = np.where(np.logical_and(df_DESeq2_all_genes[annotations]==True,df_DESeq2_all_genes['padj']<=alpha_val),"darkblue","lightgrey")
        df_DESeq2_all_genes.sort_values(by=annotations + '_color',ascending=False,inplace=True)
        ax = df_DESeq2_all_genes.plot(kind= 'scatter', x='log2FoldChange', y='negative_log10_padj', c=annotations+'_color',edgecolor='none', s=7, title = CT_i.replace('_',' '))

    ax.plot([df_DESeq2_all_genes['log2FoldChange'].min(), df_DESeq2_all_genes['log2FoldChange'].max()], [-np.log10(0.05),-np.log10(0.05)],color='black',linestyle='dotted')
    ax.plot([df_DESeq2_all_genes['log2FoldChange'].min(), df_DESeq2_all_genes['log2FoldChange'].max()], [-np.log10(0.3),-np.log10(0.3)],color='grey',linestyle='dotted')
    ax.set_xlabel('log2(FC)')
    ax.set_ylabel('-log10(adjusted p-value)')
    _,u = ax.get_ylim()
    ax.set_ylim(top=max(u,(-np.log10(alpha_val)+0.05)))
    
    if annotations=="very_small_pvals" or annotations=="small_pvals":
        if annotations=="very_small_pvals":
            th = 0.017
        else:
            th = 0.05
        if (('darkred' in np.unique(df_DESeq2_all_genes['Color_sign'].tolist())) or ('teal' in np.unique(df_DESeq2_all_genes['Color_sign'].tolist()))):
            for k, v in df_DESeq2_all_genes[['log2FoldChange','negative_log10_padj']].iterrows():
                if v['negative_log10_padj']>-np.log10(th):
                    ax.annotate(k, v, fontsize=6)
    else:
        df_sel = df_DESeq2_all_genes[df_DESeq2_all_genes[annotations]==True].copy()
        for k, v in df_sel[['log2FoldChange','negative_log10_padj']].iterrows():
            if v['negative_log10_padj']>-np.log10(0.05):
                #annotate everything true in column annotations
                ax.annotate(k, v, fontsize=8)
    ax = remove_frame_keep_axes(ax)
    if opt_save:
        os.makedirs(path_results+'figures/volcano/', exist_ok=True)
        plt.savefig(path_results+'figures/volcano/Volcano_plot_'+CT_i+'_'+str(alpha_val).replace('.','')+'.pdf', bbox_inches='tight')
    plt.close()
    plt.clf() 

def prepare_DF_for_plotting(DF_sorted_filtered,opt_aggregate_geneset_pvals_by_averaging,index_cols,p_val_cutoff):
    CT_columns = [c for c in DF_sorted_filtered.columns if c not in ["term","parentTerm"]+index_cols]
    if opt_aggregate_geneset_pvals_by_averaging:
        # macht das was es soll? unstimmigkeiten wenn man auf detailierte plots schaut!
        #make pivot table:
        DF_red_to_two_indexes = pd.pivot_table(data=DF_sorted_filtered,columns=index_cols,values =CT_columns,dropna=False,aggfunc = 'mean')#aggregates with mean but ignores nan
        # remove row if for all columns nan:
        num_cols = [x for x in DF_red_to_two_indexes.columns if x not in index_cols]
        DF_red_to_two_indexes_no_nan = DF_red_to_two_indexes.dropna(subset=num_cols,how='all')
        # remove row if none of the values is below p-val cutoff --> removes CTs
        DF_prepared = DF_red_to_two_indexes_no_nan.loc[DF_red_to_two_indexes_no_nan[num_cols].min(numeric_only=True,axis=1)<=p_val_cutoff,:]
        num_cols=[]
    else:
        #make pivot table:
        DF_sorted_p = pd.pivot_table(data=DF_sorted_filtered,index=index_cols,dropna=False)
        # remove row if for all columns nan:
        num_cols = [x for x in DF_sorted_p.columns if x not in index_cols]
        DF_sorted_p = DF_sorted_p.dropna(subset=num_cols,how='all')
        # remove row if none of the values is below p-val cutoff
        DF_prepared = DF_sorted_p.loc[DF_sorted_p[num_cols].min(numeric_only=True,axis=1)<=p_val_cutoff,:]

    return DF_prepared, num_cols

def plot_triangular_heatmaps(DF_sorted_filtered,path_info, options,n_cl,bg_str):
    if options["opt_plot_rrvgo"]:  
        GO_term_str = ['BP','MF','CC']
        for go_str in GO_term_str:
            #get GO term clustering:
            rrvgo_info_strings = ['all','up','down']
            for rrvgo_info_str in rrvgo_info_strings:
                rrvgo_filename = 'module_list_rrvgo_all_'+go_str+'_size_'+rrvgo_info_str+'.csv'
                if os.path.isfile(path_info["path_results_subdirectory"]+rrvgo_filename): 
                    module_info = pd.read_csv(path_info["path_results_subdirectory"]+rrvgo_filename)
                    plot_triangular_heatmap_rrvgo_modules(DF_sorted_filtered, module_info, go_str,options["opt_plot_rrvgo_modes"],options["p_val_cutoff"],options["opt_save"],path_info["path_results_subdirectory"], options["opt_proteomics"],options["opt_proteomics_up_and_down_sep"],n_cl, path_info["path_colors"], bg_str,rrvgo_info_str, options["opt_longread"])

    #only index group and subplot should be plotted:
    if options["opt_aggregate_geneset_pvals_by_averaging"]:
        DF_red_to_two_indexes_no_nan,_ = prepare_DF_for_plotting(DF_sorted_filtered,options["opt_aggregate_geneset_pvals_by_averaging"],options["index_cols"],options["p_val_cutoff"])
        plot_triangular_heatmap_GSEA_pathways(DF_red_to_two_indexes_no_nan, "genelists_pat_average_over_genesets", options["opt_save"], path_info["path_results_subdirectory"], options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],n_cl, path_info["path_colors"], title_str ='GSA ',background_str = bg_str)
    else:
        DF_sorted_p, num_cols = prepare_DF_for_plotting(DF_sorted_filtered,options["opt_aggregate_geneset_pvals_by_averaging"],options["index_cols"],options["p_val_cutoff"])
        for i_group,gr in enumerate(np.unique(DF_sorted_p.index.get_level_values('group'))):
            DF_tmp = DF_sorted_p.xs(gr,level = 'group')
            #TO DO: make sure all CTs are plotted no matter if all nan
            DF_tmp.dropna(subset=num_cols,how='all')
            #make sure all CTs are in heatmap always:

            # sort pathways (rows) by sign signal per ct
            if np.shape(DF_tmp)[0]>1 and gr!='synaptic-gene-ontology':
                df_tmp = DF_tmp
                df_tmp[df_tmp.isna()]=1
                pathway_list_ordered = clustermap_pathway_order(df_tmp)
                DF_tmp = DF_tmp.reindex(pathway_list_ordered)
            if DF_tmp.empty==False:
                if gr=='synaptic-gene-ontology':
                    n_pws = 60 # higher max number of pathways to make them fit on one page
                else:
                    n_pws = 35
                #if too long break down in sets of 35 pathways:
                if np.shape(DF_tmp)[0]>n_pws:
                    #print('size of '+gs_source+': '+str(np.shape(DF_p_no_nan)[0]))
                    n_rounds = int(np.ceil(np.shape(DF_tmp)[0]/n_pws))
                    #print('number of rounds: '+str(n_rounds))
                    for i_rounds in range(0,n_rounds):
                        if i_rounds<n_rounds-1:
                            DF_tmp_i = DF_tmp.iloc[i_rounds*n_pws:(i_rounds+1)*n_pws, :]
                        else:
                            DF_tmp_i = DF_tmp.iloc[i_rounds*n_pws:, :]
                        if np.shape(DF_tmp_i)[0]!=0: 
                            plot_triangular_heatmap_GSEA_pathways(DF_tmp_i, "genelists_pat_"+gr+'_'+str(i_rounds), options["opt_save"], path_info["path_results_subdirectory"],options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],options["opt_longread"],n_cl, path_info["path_colors"], title_str ='GSA ',background_str = bg_str)
                else:
                    plot_triangular_heatmap_GSEA_pathways(DF_tmp, "genelists_pat_"+gr, options["opt_save"], path_info["path_results_subdirectory"], options["opt_proteomics"], options["opt_proteomics_up_and_down_sep"],options["opt_longread"],n_cl, path_info["path_colors"], title_str ='GSA ',background_str = bg_str)

def plot_triangular_heatmap_GSEA_pathways(DF_p,geneset_name,opt_save, results_path, opt_proteomics, opt_proteomics_up_and_down_sep, opt_longread, n_clusters, path_filtered_data, title_str = '', background_str = ''):
    
    #values = create_demo_data(M, N)

    # new_cmap_greens = 'Greens'
    # new_cmap_reds = 'Reds'
    new_cmap_reds = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('darkred'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist())
    new_cmap_greens = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('teal'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                            np.array(dec_to_rgb(colors.to_rgb('white'))).tolist())
    
    if opt_proteomics:
        new_cmap_blue = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('darkgreen'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist())
    if opt_longread:
        new_cmap_blue = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('black'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                                np.array(dec_to_rgb(colors.to_rgb('white'))).tolist())
    
    #To do: 2 subplots with shared y-axis in case proteomics is up- and down-regulated together --> inly 1 column for proteomics with separate colorbar (Grays)
    cmaps = [new_cmap_reds,new_cmap_greens]  #"up","down" customize later
    #norms = [plt.Normalize(-0.5, 1) for _ in range(4)]

    max_len_label_0,max_len_label_1,max_label_length,width,add_offset_to_width,n_rows = get_settings_triang_heatmap(DF_p,results_path,opt_longread)

    if n_rows==1:
        fig_size_1, fig_size_2 = add_offset_to_width+width+2,2
    else:
        fig_size_1, fig_size_2 = add_offset_to_width+width,4+int(np.round((n_rows/7)))
    
    if (opt_proteomics == False) or (opt_proteomics == True and opt_proteomics_up_and_down_sep):

        fig, ax = plt.subplots(figsize=(fig_size_1,fig_size_2)) #,constrained_layout=True
        #fig, ax = plt.subplots(figsize=(22, 17))
        values, values_masked, values_masked_sign = get_data_from_df(DF_p,'all',n_clusters, path_filtered_data)
        #
        index_transcriptomics_1 = values[0].index
        index_transcriptomics_2 = values[1].index
        if all(index_transcriptomics_1==index_transcriptomics_2)==False:
            print('Transcriptomics indexes disagree between values[0] and values[1]!')
            return 
        #scale fontsize
        if np.shape(values)[1]>50:
            fs = 7
            fs_title = 13
        elif (np.shape(values)[1]>20 and np.shape(values)[1]<=50):
            fs = 9
            fs_title = 15
        elif opt_longread:
            fs=10
            fs_title=12
        else:
            fs = 17
            fs_title = 21
        #cmaps = [new_cmap_greens,new_cmap_reds]  # customize later
        M, N = int(np.shape(values[0])[1]), int(np.shape(values[0])[0])  # e.g. 5 columns, 4 rows
        triangul = triangulation_for_triheatmap(M, N)

        #triangle plot
        imgs = [ax.tripcolor(t, np.ravel(val), cmap=cmap, edgecolor='none',vmin=0, vmax=1)
                for t, val, cmap in zip(triangul, values_masked, cmaps)]
        #add dot inside sign triangles:
        xs, ys, _, _ = get_anchor_points_for_triangulation(M, N)
        bool_first=True
        for t,val in zip(triangul,values_masked_sign):
            #for i,id in enumerate(t.triangles):
            #    print(i.sum()/3.)
            val_vec = np.ravel(val)
            #for id,i in enumerate(val_vec):
            #    if not np.isnan(i):
            #        print(xs.ravel()[t.triangles[id]].sum()/3.)
            #        print(ys.ravel()[t.triangles[id]].sum()/3.)
            if bool_first:
                centroidX = [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                centroidY = [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                bool_first=False
            else:
                centroidX = centroidX + [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                centroidY = centroidY + [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
        plt.scatter(centroidX,centroidY, color = 'white',marker="*",s=5)

        #nan's:
        ax.patch.set(facecolor='lightgrey')
        
        ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
        ax.invert_yaxis()
        ax.margins(x=0, y=0)
        
        ##yticks should be from 0.5 to 15.5:
        #ax.set_yticks(np.arange(0, N, 1))
        #ax.set_yticks(np.arange(0.5, N+0.5, 1))
        #ax.set_yticklabels(DF_p.index,fontsize=fs, verticalalignment = 'top', horizontalalignment = 'right')
        
        labels = ['' for item in ax.get_yticklabels()]
        ax.set_yticklabels(labels)
        ax.set_ylabel('')
        #lines should be between 0 and 16
        if geneset_name.startswith("genelists_pat"):
            label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=True)
        elif geneset_name=='highlights':
            label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=False)
        elif '_rrvgo' in geneset_name:
            label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=True, add_offset=0.6)
        elif opt_longread:
            if "proteomics" in results_path:
                label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=True, add_offset=2.3, max_label_length=np.max([max_len_label_0,max_len_label_1]))
            else:
                label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=True, add_offset=0.6)
        elif geneset_name == "mitochondrial_and_ATP_pathways":
            label_group_bar_table(ax, values[0], fs, opt_add_line=True, opt_swap_indexes=True,add_offset=2.4, max_label_length=max_label_length)
        else:
            label_group_bar_table(ax, values[0], fs, opt_add_line=False, opt_swap_indexes=True,add_offset=0.2)
        
        #TO DO troubleshoot here!!
        # ax.set_xticks(np.arange(0.5, M+0.5, 1))
        ax.set_xticks(np.arange(0, M, 1))
        xticklabs = get_unique_values_in_list_while_preserving_order(values[0].columns.str.replace('based_','').str.replace(' up','').str.replace(' down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','').str.replace("_"," "))
        #_, idx = np.unique(xticklabs, return_index=True)
        # ae.set_xticklabels(xticklabs[np.sort(idx)],rotation=90, fontsize=fs, ha='right', rotation_mode='anchor')
        ax.set_xticklabels(xticklabs,rotation=90, fontsize=fs, rotation_mode='anchor', verticalalignment = 'top', horizontalalignment = 'right')
        
        #colorbars:
        fmt = '%1.2f'
        cbar_up = fig.colorbar(imgs[0], shrink =0.2, ticks=np.arange(0,0.3,0.05), pad=-0.05, format=fmt)
        cbar_down = fig.colorbar(imgs[1], shrink =0.2, ticks=np.arange(0,0.3,0.05), format=fmt)
        cbar_up.ax.tick_params(labelsize=fs)
        cbar_down.ax.tick_params(labelsize=fs)
        cbar_up.ax.set_title('up',fontsize=fs)
        cbar_down.ax.set_title('down',fontsize=fs)
        
        ax.set_aspect('equal', 'box')  # square cells
        
        # #add line before proteomics:
        if opt_proteomics:
            x = np.array([M+1, M+1])
            y = np.array([0, N])
            ax.plot(x, y,'k',linewidth = 1)
        
        # # #add lines between CT classes:
        # last_excitatory_label_id = np.max(np.where(xticklabs.str.startswith('Exc')))
        # last_inhibitory_label_id = np.max(np.where(xticklabs.str.startswith('Inh')))
        # for M in [last_excitatory_label_id,last_inhibitory_label_id]:
        #     x = np.array([M-1, M-1])
        #     y = np.array([0, N])
        #     ax.plot(x, y,'k',linewidth = 0.75)
        
        # #add lines between CT classes: 
        last_excitatory_label_id = np.max([ idx for idx,l in enumerate(xticklabs) if l.startswith('Exc')])   
        last_inhibitory_label_id = np.max([ idx for idx,l in enumerate(xticklabs) if l.startswith('Inh')])
        #last_excitatory_label_id = np.max([ idx for idx,l in enumerate(xticklabs.tolist()) if l.startswith('Exc')])
        #last_inhibitory_label_id = np.max([ idx for idx,l in enumerate(xticklabs.tolist()) if l.startswith('Inh')])
        for M in [last_excitatory_label_id,last_inhibitory_label_id]:
            x = np.array([M+1, M+1])
            y = np.array([0, N])
            ax.plot(x, y,'k',linewidth = 0.75)
            
    else:
        fig, ax = plt.subplots(1,2,gridspec_kw={'width_ratios':[int((np.shape(DF_p)[1]-3)/2)+3,1],'wspace':0, 'hspace':0},figsize=(fig_size_1,fig_size_2), sharey='row') #,constrained_layout=True
        n_sp=2
        #print(int((np.shape(DF_p)[1]-3)/2))
    
        for sp in range(0,n_sp):
            if sp==0:
                mode='transcriptomics_only'
            else:
                mode='proteomics_only'
            values,values_masked, values_masked_sign = get_data_from_df(DF_p,mode,n_clusters, path_filtered_data)
            
            if sp==0:
                index_transcriptomics_1 = values[0].index
                index_transcriptomics_2 = values[1].index
                if all(index_transcriptomics_1==index_transcriptomics_2)==False:
                    print('Transcriptomics indexes disagree between values[0] and values[1]!')
                    return 
            else:
                index_proteomics_1 = values[0].index
                index_proteomics_2 = values[1].index
                if all(index_proteomics_1==index_proteomics_2)==False:
                    print('Proteomics indexes disagree between values[0] and values[1]!')
                    return 
                if all(index_proteomics_1==index_transcriptomics_1)==False:
                    print('Proteomics indexes disagree with transcriptomics indxes!')
                    return 
            #check order of columns too:
            cols_short_1 = values[0].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')
            cols_short_2 = values[1].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')
            
            #this should never happen:
            if len(cols_short_1) < len(cols_short_2):
                #add missing columns to values[0], filled with nans for all pathways
                missing_cols = [col for col in cols_short_2 if col not in cols_short_1]
                for mc in missing_cols:
                    values[0][mc+' up'] = np.nan
                #apply same order as in values[1]
                values[0] = values[0][cols_short_2+' up']
                cols_short_1 = values[0].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')
                cols_short_2 = values[1].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')

            elif len(cols_short_2) < len(cols_short_1):
                #add missing columns to values[0], filled with nans for all pathways
                missing_cols = [col for col in cols_short_1 if col not in cols_short_2]
                for mc in missing_cols:
                    values[1][mc+' down'] = np.nan
                #apply same order as in values[0]
                values[1] = values[1][cols_short_1+' down']
                cols_short_1 = values[0].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')
                cols_short_2 = values[1].columns.str.replace('based_','').str.replace('up','').str.replace('down','').str.replace('_sum','').str.replace('_both_x','').str.replace('_both_y','').str.replace('_all','')

            if all(cols_short_1==cols_short_2)==False:
                print('Transcriptomics columns disagree between values[0] and values[1]!')
                return
            #scale fontsize
            if np.shape(values)[1]>50:
                fs = 7
                fs_title = 13
            elif (np.shape(values)[1]>20 and np.shape(values)[1]<=50):
                fs = 9
                fs_title = 15
            else:
                fs = 17
                fs_title = 21
                
            if sp==0:
                cmaps = [new_cmap_greens,new_cmap_reds]  # customize later
                M, N = int(np.shape(values[0])[1]), int(np.shape(values[0])[0])  # e.g. 5 columns, 4 rows
                triangul = triangulation_for_triheatmap(M, N)
                #triangle plot
                imgs = [ax[sp].tripcolor(t, np.ravel(val), cmap=cmap, edgecolor='none',vmin=0, vmax=1)
                        for t, val, cmap in zip(triangul, values_masked, cmaps)]
                #add dot inside sign triangles:
                xs, ys, _, _ = get_anchor_points_for_triangulation(M, N)
                bool_first=True
                for t,val in zip(triangul,values_masked_sign):
                    #for i,id in enumerate(t.triangles):
                    #    print(i.sum()/3.)
                    val_vec = np.ravel(val)
                    #for id,i in enumerate(val_vec):
                    #    if not np.isnan(i):
                    #        print(xs.ravel()[t.triangles[id]].sum()/3.)
                    #        print(ys.ravel()[t.triangles[id]].sum()/3.)
                    if bool_first:
                        centroidX = [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        centroidY = [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        bool_first=False
                    else:
                        centroidX = centroidX + [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        centroidY = centroidY + [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                plt.scatter(centroidX,centroidY, color = 'white',marker="*",s=5)
                #nan's:
                #ax[sp].patch.set(hatch='///',edgecolor='black',facecolor='ivory')
                #nan's:
                ax[sp].patch.set(facecolor='lightgrey')
                
                #does this revert everything, also tick labels??
                #ax[sp].invert_yaxis()
                
                ax[sp].grid(which='major', color='#DDDDDD', linewidth=0.8)
                ax[sp].margins(x=0, y=0)
                
                labels = ['' for item in ax[sp].get_yticklabels()]
                ax[sp].set_yticklabels(labels)
                ax[sp].set_ylabel('')
                #lines should be between 0 and 16
                
                if geneset_name.startswith("genelists_pat"):
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=True, opt_swap_indexes=True,add_offset=0.5,max_label_length=max_label_length)
                # elif geneset_name.startswith("mitoxplorer"):
                #     label_group_bar_table(ax[sp], DF_p, fs, opt_add_line=False, opt_swap_indexes=False,add_offset=0.3)
                elif '_rrvgo' in geneset_name:
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=True, opt_swap_indexes=True,add_offset=0.6,max_label_length=max_label_length)
                elif geneset_name.startswith("GO_BP"):
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=False, opt_swap_indexes=False,add_offset=0.4,max_label_length=max_label_length)
                elif geneset_name=='highlights':
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=True, opt_swap_indexes=True,add_offset=0.5,max_label_length=90)
                elif geneset_name=="mitochondrial_and_ATP_pathways":
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=True, opt_swap_indexes=True,add_offset=0.4,max_label_length=max_label_length)
                else:
                    label_group_bar_table(ax[sp], values[0], fs, opt_add_line=False, opt_swap_indexes=False,add_offset=0.4,max_label_length=max_label_length)
                
                
                
                ax[sp].set_xticks(np.arange(0, M, 1))
                ax[sp].grid(which='major', color='#DDDDDD', linewidth=0.8)
                # ax[sp].set_yticks(np.arange(0, N, 1))
                # ax[sp].set_yticklabels(DF_p.index,fontsize=fs, verticalalignment = 'top', horizontalalignment = 'right')
                
                # ax.set_yticks(np.arange(0.5, N+0.5, 1))
                labels = ['' for item in ax[sp].get_yticklabels()]
                ax[sp].set_yticklabels(labels)
                ax[sp].set_ylabel('')
                
                ax[sp].invert_yaxis()
                    
                #locs = ax[sp].get_yticks()
                
                #values for colorbar:
                imgs_up = imgs[0]
                imgs_down = imgs[1]
                # M_first_sp = M
            else:
                cmaps = [new_cmap_blue,new_cmap_blue]  # customize later
                M = int(np.shape(values[0])[1])  #update M
                
                triangul = triangulation_for_triheatmap(M, N)
                #heatmap
                imgs = [ax[sp].tripcolor(t, np.ravel(val), cmap=cmap, edgecolor='none',vmin=0, vmax=1)
                        for t, val, cmap in zip(triangul, values_masked, cmaps)]
                
                #add dot inside sign triangles:
                xs, ys, _, _ = get_anchor_points_for_triangulation(M, N)
                bool_first=True
                for t,val in zip(triangul,values_masked_sign):
                    #for i,id in enumerate(t.triangles):
                    #    print(i.sum()/3.)
                    val_vec = np.ravel(val)
                    #for id,i in enumerate(val_vec):
                    #    if not np.isnan(i):
                    #        print(xs.ravel()[t.triangles[id]].sum()/3.)
                    #        print(ys.ravel()[t.triangles[id]].sum()/3.)
                    if bool_first:
                        centroidX = [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        centroidY = [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        bool_first=False
                    else:
                        centroidX = centroidX + [xs.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                        centroidY = centroidY + [ys.ravel()[t.triangles[id]].sum()/3.  for id,i in enumerate(val_vec) if not np.isnan(i)]
                plt.scatter(centroidX,centroidY, color = 'white',marker="*",s=5)

                #nan's:
                #ax[sp].patch.set(hatch='///',edgecolor='black',facecolor='ivory')
                ax[sp].patch.set(facecolor='lightgrey')
                #imgs_prot = ae.imshow(values[0],cmap=new_cmap_grays)
                ax[sp].set_xticks(np.arange(0, 1, 1))
                #set the same yticks as in first subplot
                #ax[sp].set_yticks(locs)
                #ax[sp].set_yticklabels([])

                # M_second_sp = M
                #ax[sp].set_aspect('equal', 'box')  # square cells

            #_, idx = np.unique(xticklabs, return_index=True)
            # ae.set_xticklabels(xticklabs[np.sort(idx)],rotation=90, fontsize=fs, ha='right', rotation_mode='anchor')
            ax[sp].set_xticklabels(cols_short_1,rotation=90, fontsize=fs, rotation_mode='anchor', verticalalignment = 'top', horizontalalignment = 'right')   
        
            if sp==0:
                # #add lines between CT classes:
                last_excitatory_label_id = np.max([ idx for idx,l in enumerate(cols_short_1.tolist()) if l.startswith('Exc')])
                last_inhibitory_label_id = np.max([ idx for idx,l in enumerate(cols_short_1.tolist()) if l.startswith('Inh')])
                for M in [last_excitatory_label_id,last_inhibitory_label_id]:
                    x = np.array([M+1, M+1])
                    y = np.array([0, N])
                    ax[0].plot(x, y,'k',linewidth = 0.75)
        
    if len(title_str)>0:
        add_str = title_str.replace(' ','_')
        fig.suptitle(title_str+background_str+' result of '+geneset_name, fontsize = fs_title)
    else:
        add_str = 'GSEA_'
        fig.suptitle('GSEA result of '+geneset_name, fontsize = fs_title)
    #fig.tight_layout()
    #plt.tight_layout()
    #plt.show() 
    
    if opt_proteomics==True:
        proteomics_str = '_with_proteomics'
        if opt_proteomics_up_and_down_sep==False:
            # make cells quadratic:
            #for sp in range(0,n_sp):
            #    ax[sp].set_aspect('equal', 'box')  # square cells
            plt.margins(x=0, y=0, tight=True)
            
            #colorbars:
            fmt = '%1.2f'
            
            # Create new axes according to image position
            cax_up = fig.add_axes([ ax[1].get_position().x0+0.05,
                                    ax[1].get_position().y0,
                                    0.02,
                                    ax[1].get_position().height-0.05])
            cax_down = fig.add_axes([ ax[1].get_position().x0+0.15,
                                    ax[1].get_position().y0,
                                    0.02,
                                    ax[1].get_position().height-0.05])
            cax_prot = fig.add_axes([ ax[1].get_position().x0+0.25,
                                    ax[1].get_position().y0,
                                    0.02,
                                    ax[1].get_position().height-0.05])
    
            cbar_up = plt.colorbar(imgs_up, shrink =0.2, ticks=np.arange(0,0.3,0.05), pad=-0.5, format=fmt, cax=cax_up) #cax=cax1)
            cbar_down = plt.colorbar(imgs_down, shrink =0.2, ticks=np.arange(0,0.3,0.05), format=fmt, cax=cax_down) #cax=cax2)
            cbar_prot = plt.colorbar(imgs[0], shrink =0.2, ticks=np.arange(0,0.3,0.05), format=fmt, cax=cax_prot) #cax=cax3)
            cbar_up.ax.tick_params(labelsize=fs)
            cbar_down.ax.tick_params(labelsize=fs)
            cbar_prot.ax.tick_params(labelsize=fs)
            cbar_up.ax.set_title('up',fontsize=fs)
            cbar_down.ax.set_title('down',fontsize=fs)
            cbar_prot.ax.set_title('proteomics',fontsize=fs)
        else:
            proteomics_str = proteomics_str + '_up_down_sep'
    else:
        plt.margins(x=0, y=0, tight=True)
        proteomics_str = ''
    #plt.show()   
    
    #more fine tuning for legend: https://stackoverflow.com/questions/43697234/making-colorbar-values-integer-in-a-heatmap-matplotlib-seaborn
    print(results_path)
    if opt_save:
        path=results_path+'figures/heatmaps/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        fig.savefig(path + 'HM3_deregGenesP_'+add_str+geneset_name.replace(':','_')+proteomics_str+'.pdf',bbox_inches='tight')
        plt.close(fig)
        plt.clf() 

def get_settings_triang_heatmap(DF_p,results_path,opt_longread):
    n_rows = np.shape(DF_p)[0]
    if "Mitochondria" in results_path:
        max_len_label_0 = 100
        max_len_label_1 = 100
    else:
        max_len_label_0 = np.min([100,np.max(DF_p.index.get_level_values(0).str.len())])
        if "proteomics" in results_path or DF_p.index.nlevels: #only one index
            max_len_label_1 = 100
        else:
            max_len_label_1 = np.min([100,np.max(DF_p.index.get_level_values(1).str.len())])
    if DF_p.index.nlevels == 1:
        width = 4+int(np.round((max_len_label_0)/5))
    elif DF_p.index.nlevels == 2:
        width = 4+int(np.round((max_len_label_0+max_len_label_1)/10))
    else:
        print('This is not implemented for more than 2 indices yet! Go to plot_triangular_heatmap_GSEA_pathways() in self_fun_plot.py!')
    add_offset_to_width=0
    if n_rows<20:
        add_offset_to_width += 2
        if n_rows <=16:
            add_offset_to_width += 1
            if n_rows <10:
                add_offset_to_width += 1
    if opt_longread:
        add_offset_to_width += 5
    if "Mitochondria" in results_path:
        max_label_length = 100
    else:
        max_label_length = 100

    return max_len_label_0,max_len_label_1,max_label_length,width,add_offset_to_width,n_rows

def get_anchor_points_for_triangulation(M, N):
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)
    trianglesNW = [(i + j*(M+1), i+1 + j*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    trianglesSE = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]

    return xs, ys, trianglesNW, trianglesSE

def triangulation_for_triheatmap(M, N):
    xs, ys, trianglesNW, trianglesSE = get_anchor_points_for_triangulation(M, N)

    return [Triangulation(xs.ravel(), ys.ravel(), triangles) for triangles in [trianglesNW,trianglesSE]]

def plot_triangular_heatmap_rrvgo_modules(DF_sorted_filtered, module_info, go_str,opt_plot_rrvgo_modes,p_val_cutoff,opt_save,folder, opt_proteomics,opt_proteomics_up_and_down_sep,n_cl, path_colors, bg_str, rrvgo_info_str, opt_longread):

    #for DF_sorted_filtered get selection of current go term:
    if 'subgroup' in DF_sorted_filtered.index.names:
        DF_sorted_filtered_sel=DF_sorted_filtered[DF_sorted_filtered.index.get_level_values('subgroup')=='GO:'+go_str]    
    else:
        DF_sorted_filtered_sel=DF_sorted_filtered[DF_sorted_filtered['subgroup']=='GO:'+go_str]
    
    #split geneset in go id and pathway name
    DF_sorted_filtered_sel.reset_index(inplace=True)
    DF_sorted_filtered_sel[['go','pw']] = DF_sorted_filtered_sel['geneset'].str.rsplit(pat=r' ',n=1,expand=True)
        
    #merge dataframes:
    DF_modules = pd.merge(DF_sorted_filtered_sel,module_info,on='go')
    
    #get rid of unnecessary columns:
    DF_modules.drop(columns=['group', 'subgroup', 'geneset','go', 'pw', 'Unnamed: 0', 'cluster','parent', 'parentSimScore', 'score', 'size'],inplace=True)  
    
    
    #To Do: make sure every parent appears in plot/DF
    for opt_plot_rrvgo_mode in opt_plot_rrvgo_modes:
        if opt_plot_rrvgo_mode=='average_for_each_parent':
            average_str = '_averaged'
            opt_aggregate_geneset_pvals_by_averaging=True
            DF_modules_red_to_two_indexes_no_nan,_ = prepare_DF_for_plotting(DF_modules,opt_aggregate_geneset_pvals_by_averaging,['parentTerm'],p_val_cutoff)
            opt_proteomics = False
            #transpose data frame:
            DF_modules_red_to_two_indexes_no_nan = DF_modules_red_to_two_indexes_no_nan.transpose()
        elif opt_plot_rrvgo_mode=='significant_sorted_by_parent':
            average_str = ''
            opt_aggregate_geneset_pvals_by_averaging=False
            DF_modules_red_to_two_indexes_no_nan,_ = prepare_DF_for_plotting(DF_modules,opt_aggregate_geneset_pvals_by_averaging,['parentTerm', 'term'],p_val_cutoff)
       
        

        #if too long break down in sets of n_pw_max pathways:
        n_pw_max = 45
        if np.shape(DF_modules_red_to_two_indexes_no_nan)[0]>n_pw_max:
            #print('size of '+gs_source+': '+str(np.shape(DF_p_no_nan)[0]))
            n_rounds = int(np.ceil(np.shape(DF_modules_red_to_two_indexes_no_nan)[0]/n_pw_max))
            #print('number of rounds: '+str(n_rounds))
            for i_rounds in range(0,n_rounds):
                if i_rounds<n_rounds-1:
                    DF_modules_red_to_two_indexes_no_nan_i = DF_modules_red_to_two_indexes_no_nan.iloc[i_rounds*n_pw_max:(i_rounds+1)*n_pw_max, :]
                else:
                    DF_modules_red_to_two_indexes_no_nan_i = DF_modules_red_to_two_indexes_no_nan.iloc[i_rounds*n_pw_max:, :]
                if np.shape(DF_modules_red_to_two_indexes_no_nan_i)[0]!=0: 
                    plot_triangular_heatmap_GSEA_pathways(DF_modules_red_to_two_indexes_no_nan_i, "GO_"+go_str+'_rrvgo_'+rrvgo_info_str+'_'+average_str+'_'+str(i_rounds), opt_save, folder, opt_proteomics,opt_proteomics_up_and_down_sep,opt_longread, n_cl, path_colors,title_str ='GSA ',background_str = bg_str)
        else:
            plot_triangular_heatmap_GSEA_pathways(DF_modules_red_to_two_indexes_no_nan, "GO_"+go_str+'_rrvgo_'+rrvgo_info_str+'_'+average_str, opt_save, folder, opt_proteomics,opt_proteomics_up_and_down_sep,opt_longread,n_cl, path_colors,title_str ='GSA ',background_str = bg_str)

        #plot triangle plot for 'parentTerm', 'term'


def remove_frame_keep_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    return ax


def loadPltSettings(fontSize,markerSize):
    plt.style.use('classic')
    plt.box(False)
    plt.rcParams['grid.alpha'] = 0
    plt.rcParams['font.size'] = fontSize
    plt.rcParams['axes.edgecolor'] ='black'
    plt.rcParams['axes.facecolor'] ='black'#'white'
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['figure.edgecolor'] = 'white'
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['figure.titlesize'] = fontSize
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.labelsize']= fontSize
    plt.rcParams['legend.fontsize'] = fontSize
    plt.rcParams['boxplot.flierprops.markersize']= markerSize
    plt.rcParams["axes.formatter.use_mathtext"]=True
    plt.rcParams["axes.formatter.limits"] = (0,0)
    plt.rcParams["legend.numpoints"] = 1
    plt.rcParams["legend.scatterpoints"] = 1
    plt.rcParams["figure.figsize"] = [11,8]
    #remove frame
    #plt.rcParams['axes.spines.right'] = False
    #plt.rcParams['axes.spines.top'] = False
    return plt

### calculate functions:
def calculate_percentage_random_DEGs(df_true,path_results_short,design_folder_str,cell_type_folder_str,CT_i,alpha_val):
    df_true['log2FoldChange_abs'] = df_true['log2FoldChange'].abs()
    df_true.sort_values('log2FoldChange_abs',axis=0,ascending=False,inplace=True)
    #remove whitespaces
    df_true['Gene'] = df_true['Gene'].str.strip()
    df_true.set_index('Gene',inplace=True)
    df_true['P_random_DEG'] = 0
    #use order of genes to sort percentages of permuted datasets where gene_i is wrongly identified as DEG:
    
    #load results from randomized group label data:
    path_randomized_results = path_results_short+'randomized_group_labels/'+design_folder_str+str.split(cell_type_folder_str,'/')[-1]+'/'
    cp = os.getcwd()
    bool_first_df_rand = True
    if os.path.exists(path_randomized_results):
        os.chdir(path_randomized_results)
        for i in range(1,101):
            filename = 'df_results_no_shrinkage'+str(i)+'.csv'
            if os.path.exists(filename):
                df_rand = pd.read_csv(filename)
                df_rand['Unnamed: 0'] = df_rand['Unnamed: 0'].str.strip()
                df_rand.set_index('Unnamed: 0',inplace=True)
                #resort according to df_true gene order:
                df_rand = df_rand.reindex(index=df_true.index.tolist())
                if np.shape(df_rand[df_rand.padj<=alpha_val])[0]>0:
                    #check which genes fullfill padj < alpha_val
                    df_true.loc[(df_rand.padj<=alpha_val),'P_random_DEG'] += 1 
                    df_rand.drop(df_rand[(df_rand.padj > alpha_val) | (df_rand.padj.isnull())].index, inplace=True)
                    #put in 1 large df:
                    if bool_first_df_rand:
                        df_rand_sign = df_rand
                        bool_first_df_rand = False
                    else:
                        df_rand_sign = pd.concat([df_rand_sign,df_rand])
            #add ct info to df_true
            df_true['Celltype'] = CT_i
    
    #not interested in genes that were not significant in true labeling or had nan adj p-val: 
    df_true.drop(df_true[(df_true.padj > alpha_val) | (df_true.padj.isnull())].index, inplace=True)
    #might result in empty data_frame
    if bool_first_df_rand==False:
        df_rand_sign['Celltype'] = CT_i
    else:
        df_rand_sign = pd.DataFrame([], columns=['nothing'])
    #go back to initial path:
    os.chdir(cp)
    return df_true, df_rand_sign


def dec_to_rgb(rgb_tuple_dec):
    r = rgb_tuple_dec[0]*255.0
    b = rgb_tuple_dec[1]*255.0
    g = rgb_tuple_dec[2]*255.0
    return (r,b,g)
    #https://stackoverflow.com/questions/25404998/how-can-i-convert-an-rgb-input-into-a-decimal-color-code
    #return (r << 16) + (g << 8) + b