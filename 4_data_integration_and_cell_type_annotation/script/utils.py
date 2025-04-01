# author: Lisa Bast
# date: Created on Mon Nov 16 15:14:43 2020
# version: 0.0.1
# about: utils for cell type clustering 
#        create csv files from downloaded reference brain map to use it with scMatch automated cell type annotation

import os
import pandas as pd
import anndata 
import loompy
import numpy as np
import scanpy as sc
import gc
import random
import plotly.graph_objects as go
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics 
from sklearn.decomposition import PCA

def get_paths(path_project,opt_output):
    
    if opt_output == "Cell type annotation":
        path_data = path_project + "3_quality_control/output/filtered/"
        path_code = path_project + "4_data_integration_and_cell_type_annotation/script/"    
        path_results = path_project + "4_data_integration_and_cell_type_annotation/output/CT_anno_scmap/"
    elif opt_output == "Cell type-specific files":
        path_data = path_project + "4_data_integration_and_cell_type_annotation/output/CT_anno_scmap/"
        path_code = path_project + "4_data_integration_and_cell_type_annotation/script/"    
        path_results = path_project + "4_data_integration_and_cell_type_annotation/output/CT_specific_files/"
    return path_data, path_code, path_results

def get_sampleIDS_nCells_diseaseStatus_sexStatus(D):
    sample_IDs = np.unique(D.ca.Donor)
    n_cells = np.empty(len(sample_IDs), dtype=int)
    disease_status = list()
    sex_status = list()
    for d_id,donor in enumerate(sample_IDs):
        # number of cells per sample
        n_cells[d_id] = len(np.unique(D.ca.CellID[np.where(D.ca.Donor==donor)]))
        disease_status.append(np.unique(D.ca.Disease[np.where(D.ca.Donor==donor)])[0])
        sex_status.append(np.unique(D.ca.Sex[np.where(D.ca.Donor==donor)])[0])
    return sample_IDs,n_cells,disease_status,sex_status

def get_CT_assigned(col_attr,opt_celltype_groups,ct_col_of_interest):
    if opt_celltype_groups == 'scmap':
        CT_list = np.unique(col_attr[ct_col_of_interest]).tolist()
        CT_list.remove('nan')
        CT_list.remove('unassigned')
    elif opt_celltype_groups == 'conos_cluster_based':
        CT_list = np.unique(col_attr[ct_col_of_interest]).tolist()
        if 'removed' in CT_list:
            CT_list.remove('removed')
        if 'none (removed)' in CT_list:
            CT_list.remove('none (removed)')
    return CT_list

def add_cell_class_info(df,CT_col_name):
    df['cell class'] = df[CT_col_name]#['non-neuronal']*len(GT['cell type'])
    for s_id in range(0,np.shape(df)[0]):
    #for s_id,sample in enumerate(df['sample'].tolist()):
        #print(s_id)
        if not isinstance(df['cell class'][s_id],str):
            if np.isnan(df['cell class'][s_id]):
                df['cell class'][s_id] = 'NaN'
        else:
            if df['cell class'][s_id].startswith('Exc'):
                df['cell class'][s_id] = 'Exc'
            elif df['cell class'][s_id].startswith('Inh'):
                df['cell class'][s_id] = 'Inh'
            elif df['cell class'][s_id].startswith('unassigned'):
                df['cell class'][s_id] = 'unassigned'
    return df 

def get_DF_annotation_result_test(GT, red_methods, tool, n_genes_ref, path_results_files):
    first_time_bool = True
    for m,method in enumerate(red_methods):
        for n in n_genes_ref:
            file_str_cluster_score = 'R_cluster_score_'+str(n)+'_features_based_on_'+method+'.csv'
            file_str_cluster_ann = 'R_cluster_ann_CT_'+str(n)+'_features_based_on_'+method+'.csv'
            file_str_cell_score = 'R_cell_score_'+str(n)+'_features_based_on_'+method+'.csv'
            file_str_cell_ann = 'R_cell_ann_cellID_'+str(n)+'_features_based_on_'+method+'.csv'
            file_str_cell2cluster_score = 'R_cell2cluster_score_'+str(n)+'_features_based_on_'+method+'.csv'
            file_str_cell2cluster_ann = 'R_cell2cluster_ann_CT_'+str(n)+'_features_based_on_'+method+'.csv'
            os.chdir(path_results_files)
            df_ann = pd.read_csv(file_str_cluster_ann)
            df_ann = add_cell_class_info(df_ann,'ABM')
            CT_ann_1 = df_ann['ABM']
            Bools_correct_CT_1 = GT['cell type']==df_ann['ABM']
            Bools_correct_CC_1 = GT['cell class']==df_ann['cell class']
            cohens_k_1 = sklearn.metrics.cohen_kappa_score(CT_ann_1.tolist(),GT['cell type'].tolist())
            if first_time_bool:
                A_bools_str = [tool + '_cluster_' + method]
                N_genes = [n]
            else:
                A_bools_str.append(tool + '_cluster_' + method)
                N_genes.append(n)
            Bools_unassigned_1 = df_ann['ABM']=='unassigned'
            Bools_unassigned_CC_1 = df_ann['cell class']=='unassigned'
            #A_bools_str.append(tool + '_cluster_' + method + '_' +str(n)+ 'genes_UA')
            #N_genes.append(n)
            
            df_ann = pd.read_csv(file_str_cell2cluster_ann)
            df_ann = add_cell_class_info(df_ann,'ABM')
            CT_ann_2 = df_ann['ABM']
            Bools_correct_CT_2 = GT['cell type']==df_ann['ABM']
            Bools_correct_CC_2 = GT['cell class']==df_ann['cell class']
            cohens_k_2 = sklearn.metrics.cohen_kappa_score(CT_ann_2.tolist(),GT['cell type'].tolist())
            A_bools_str.append(tool + '_cell2cluster_' + method )
            N_genes.append(n)
            Bools_unassigned_2 = df_ann['ABM']=='unassigned'
            Bools_unassigned_CC_2 = df_ann['cell class']=='unassigned'
            if first_time_bool:
                A_bools_correct_CTs = np.vstack([Bools_correct_CT_1,Bools_correct_CT_2])#,Bools_correct_CT_3])
                A_bools_correct_CCs = np.vstack([Bools_correct_CC_1,Bools_correct_CC_2])#,Bools_correct_CC_3])
                A_bools_unassigned_CTs = np.vstack([Bools_unassigned_1,Bools_unassigned_2])#,Bools_unassigned_3])
                A_bools_unassigned_CCs = np.vstack([Bools_unassigned_CC_1,Bools_unassigned_CC_2])#,Bools_CC_unassigned_3])
                A_CT_ann = np.vstack([CT_ann_1,CT_ann_2])#,CT_ann_3.flatten()])
                A_cohens_k = np.vstack([cohens_k_1,cohens_k_2])#,cohens_k_3])
                first_time_bool = False
            else:
                A_bools_correct_CTs = np.vstack([A_bools_correct_CTs,Bools_correct_CT_1,Bools_correct_CT_2])#,Bools_correct_CT_3])
                A_bools_correct_CCs = np.vstack([A_bools_correct_CCs,Bools_correct_CC_1,Bools_correct_CC_2])#,Bools_correct_CC_3])
                A_bools_unassigned_CTs = np.vstack([A_bools_unassigned_CTs,Bools_unassigned_1,Bools_unassigned_2])#,Bools_unassigned_3])
                A_bools_unassigned_CCs = np.vstack([A_bools_unassigned_CCs,Bools_unassigned_CC_1,Bools_unassigned_CC_2])#,Bools_unassigned_CC_3])
                A_CT_ann = np.vstack([A_CT_ann,CT_ann_1,CT_ann_2])#,CT_ann_3.flatten()])
                A_cohens_k = np.vstack([A_cohens_k,cohens_k_1,cohens_k_2])#,cohens_k_3])
    DF_R = pd.DataFrame(data = {'bools_correct_cell_type': np.sum(A_bools_correct_CTs,axis=1),'bools_correct_cell_class': np.sum(A_bools_correct_CCs,axis=1), 'bools_unassigned_cell_type': np.sum(A_bools_unassigned_CTs,axis=1),'bools_unassigned_cell_class': np.sum(A_bools_unassigned_CCs,axis=1),'n_genes': N_genes, 'method': A_bools_str})
    DF_R['Cohens_kappa'] = A_cohens_k
    N_cells_test = np.shape(A_bools_correct_CTs)[1]
    DF_R['P_correct'] = DF_R['bools_correct_cell_type']/N_cells_test*100
    DF_R['P_unassigned'] = DF_R['bools_unassigned_cell_type']/N_cells_test*100
    DF_R['P_false'] = 100-DF_R['P_correct']-DF_R['P_unassigned']
    DF_R['P_CC_correct'] = DF_R['bools_correct_cell_class']/N_cells_test*100
    DF_R['P_CC_unassigned'] = DF_R['bools_unassigned_cell_class']/N_cells_test*100
    DF_R['P_CC_false'] = 100-DF_R['P_CC_correct']-DF_R['P_CC_unassigned']
    return DF_R
def format_ref_data(main_path, ref_data_set, ref_data_set_short, n_CT_vec, opt_all_genes, red_method, test_data_path, test_data_filename, IDs_cell_sel, samples_selected, opt_ref_data_as_query_data, n_variable_genes):
    
    if ref_data_set == 'Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq':
        if opt_ref_data_as_query_data:
            filename = 'matrix_46917_cells_GT_train.csv'
        else:
            filename = 'metadata.csv'
    else:
        print('The specified ref data set is not available!')
    #read files
    path_raw_ref_files = main_path+'/data/reference_data_sets/' +ref_data_set_short
    os.chdir(path_raw_ref_files)
    
    if ref_data_set == 'Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq' and opt_ref_data_as_query_data:
        aD = anndata.read_csv("matrix_46917_cells_train.csv")
        aD_red_cols = aD[aD.obs.group=='train',:]
    else:
        print('The specified ref data set is not available!')
   
    genes_priority = ['SLC17A7','SLC32A1','MBP','AQP4','AIF1','PDGFRA','GAD1','ADARB2','LAMP5','PAX6','VIP','LHX6','SST','PVALB','CUX2','RORB','THEMIS','FEZF2','CTGF','MOG','CLDN5','MUSTN1','FYB','OPALIN','FYB']
    if opt_all_genes:
        #delete all genes that are not present in test data set anyways to reduce file size
        os.chdir(test_data_path)
        D_test = loompy.connect(test_data_filename)
        genes_test = D_test.ra['Gene']
        D_test.close()
        if len(genes_test)<len(aD_red_cols.var_names):
            genes_intersect = np.intersect1d(aD_red_cols.var_names,genes_test)
            aD_red = aD_red_cols[:, genes_intersect]
            df_symbol = aD_red.copy().transpose().to_df()
            add_str = '_red'
        else:
            df_symbol = aD_red_cols.copy().transpose().to_df()
            add_str = ''
    else:
        #determine top k variable genes based on ref data and keep only those
        if red_method == 'highest_dropout_rate':
            os.chdir(main_path+'/data/reference_data_sets/'+ref_data_set_short+'/')
            df_HDO = pd.read_csv('list_genes_selected_high_dropout_'+str(n_variable_genes)+'.csv')
            add_str = '_red_to_'+str(n_variable_genes)+'_HDO'
            #ab hier weiter debuggen:
            genes_with_high_dropout_rate = df_HDO['x'].tolist()
            aD_red = aD_red_cols[:,np.logical_or(aD_red_cols.var_names.isin(genes_with_high_dropout_rate), aD_red_cols.var_names.isin(genes_priority))]   
        elif red_method.startswith('HVG'):
            if red_method == 'HVG_with_seurat_v3':
                sc.pp.highly_variable_genes(aD_red_cols,flavor='seurat_v3',n_top_genes=n_variable_genes,inplace=True)
                hvgs = aD_red_cols.var.highly_variable
            elif red_method == 'HVG_with_cell_ranger':
                aD_red_cols_log1p = sc.pp.log1p(aD_red_cols, copy=True)
                sc.pp.highly_variable_genes(aD_red_cols_log1p,flavor='cell_ranger',n_top_genes=n_variable_genes,inplace=True)
                hvgs = aD_red_cols_log1p.var.highly_variable
            else:
                print('red_method '+red_method+' not implemented!')
            aD_red = aD_red_cols[:,np.logical_or(hvgs, aD_red_cols.var_names.isin(genes_priority))]   
            add_str = '_red_to_'+str(n_variable_genes)+'_'+red_method                              
        else:
            print('red_method '+red_method+' not implemented!') 
        df_symbol = aD_red.copy().transpose().to_df()
    #delete all genes with less than 20 counts in total
    G_kept = df_symbol.sum(axis=1)>20
    G_kept_array = np.intersect1d(df_symbol.index.values,G_kept[G_kept].index.values)
    df_symbol = df_symbol.loc[G_kept_array.tolist()]

    #rename first column: symbol
    df_symbol = df_symbol.reset_index()
    df_symbol.rename(columns = {'index':'symbol'},inplace=True)
    path_formatted_ref_files = main_path+'data/reference_data_sets/'+ref_data_set_short+'/'+add_str[1:]
    os.chdir(path_formatted_ref_files)
    df_symbol.to_csv('9606_symbol.csv',index=False)

    #create and save map:
    for n_CTs in n_CT_vec:
        df_map = pd.read_csv(filename)
        df_map.drop(columns=['index_whole_data_set'],inplace=True)
        os.chdir(path_formatted_ref_files)
        df_map.to_csv('9606_map_'+str(n_CTs)+'.csv',index=True)
    
def merge_CT_cluster_of_metadata(ref_data_set,metadata,n_CTs):
    if ref_data_set == 'Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq':
        metadata.rename(columns = {'cell_type_alias_label':'cell_type_alias_label_120CTs'},inplace=True)
        if n_CTs!=18:
            metadata['cell_type_alias_label_76CTs'] = metadata['cell_type_alias_label_120CTs'].copy()
            old_names = [['Inh L1 LAMP5 GGT8P','Inh L1 LAMP5 NDNF'],
                            ['Inh L1-4 LAMP5 DUSP4','Inh L6 LAMP5 C1QL2'],
                            ['Inh L1-6 LAMP5 CA13', 'Inh L5-6 LAMP5 SFTA3'],
                            ['Inh L1 PAX6 CA4','Inh L1 PAX6 GRIP2'],
                            ['Inh L1-3 PAX6 NABP1','Inh L1-6 VIP RCN1'],
                            ['Inh L1 VIP PRSS8','Inh L1 VIP TNFAIP8L3'],
                            ['Inh L1 ADARB2 ADAM33','Inh L1 SST CXCL14'],
                            ['Inh L1 ADARB2 DISP2','Inh L1 VIP SOX11'],
                            ['Inh L1-6 VIP PENK','Inh L1-5 VIP KCNJ2'],
                            ['Inh L1 VIP PCDH20','Inh L1-2 VIP PPAPDC1A'],
                            ['Inh L2-6 VIP VIP','Inh L3-6 VIP KCTD13'],
                            ['Inh L1-2 VIP RPL41P3','Inh L1-3 VIP ACHE'],
                            ['Inh L2-4 VIP LGI2','Inh L1-4 VIP CHRNA2'],
                            ['Inh L2-5 VIP TOX2','Inh L2-4 VIP DSEL'],
                            ['Inh L1-3 VIP ZNF322P1','Inh L3 VIP CBLN1'],
                            ['Inh L5-6 SST ISOC1','Inh L5-6 SST KLHL14'],
                            ['Inh L3-5 SST MAFB','Inh L4-6 SST MTHFD2P6'],
                            ['Inh L6 LHX6 GLP1R','Inh L5-6 PVALB FAM150B'],
                            ['Inh L5 PVALB CNTNAP3P2','Inh L2-4 PVALB C8orf4'],
                            ['Inh L1-2 PVALB TAC1','Inh L1-3 PVALB WFDC2'],
                            ['Inh L1-6 PVALB SCUBE3','Inh L3-6 PVALB MFI2'],
                            ['Exc L2-4 RORB GRIK1','Exc L3-4 RORB RPS3P6'],
                            ['Exc L4 RORB BHLHE22','Exc L4 RORB CCDC168'],
                            ['Exc L4-5 RORB ASCL1','Exc L4-5 RORB AIM2'],
                            ['Exc L3 LINC00507 PSRC1','Exc L3-5 LINC00507 SLN'],
                            ['Exc L3 RORB CARTPT','Exc L3-4 RORB FOLH1B'],
                            ['Exc L3-4 RORB SEMA6D','Exc L3-4 RORB PRSS12'],
                            ['Exc L3 LINC00507 CTXN3','Exc L3 THEMIS PLA2G7'],
                            ['Exc L5-6 THEMIS OR1J1','Exc L6 THEMIS EGR3'],
                            ['Exc L6 THEMIS LINC00343','Exc L5-6 THEMIS TMEM233'],
                            ['Exc L3-5 RORB CMAHP','Exc L3-5 RORB CD24'],
                            ['Exc L5-6 RORB LINC00320','Exc L5 RORB LINC01202'],
                            ['Exc L4-5 RORB HNRNPA1P46','Exc L4-5 RORB RPL31P31'],
                            ['Exc L3-5 THEMIS UBE2F','Exc L5 RORB SNHG7'],
                            ['Exc L6 THEMIS C6orf48','Exc L5-6 THEMIS GPR21'],
                            ['Exc L5-6 FEZF2 ANKRD20A1','Exc L6 FEZF2 FAM95C'],
                            ['Exc L6 FEZF2 TBC1D26','Exc L6 FEZF2 KRT17'],
                            ['Exc L6 FEZF2 SLITRK6','Exc L6 FEZF2 TBCC'],
                            ['Exc L5 FEZF2 SCN7A','Exc L5 FEZF2 MORN2'],
                            ['Exc L5 FEZF2 DYRK2','Exc L5-6 FEZF2 MYBPHL'],
                            ['Exc L5-6 FEZF2 CYP26B1','Exc L5-6 FEZF2 RSAD2'],
                            ['Astro L1 FGFR3 MT1G','Astro L1 FGFR3 FOS'],
                            ['Oligo L4-6 MOBP COL18A1','Oligo L4-6 OPALIN'],
                            ['VLMC L1-3 CYP1B1','Peri L1-6 MUSTN1']]
            new_names = ['Inh L1 LAMP5 GGT8P or NDNF',
                            'Inh L1-4 LAMP5 DUSP4 or L6 LAMP5 C1QL2',
                            'Inh L1-6 LAMP5 CA13 or L5-6 LAMP5 SFTA3',
                            'Inh L1 PAX6 CA4 or GRIP2',
                            'Inh L1-3 PAX6 NABP1 or L1-6 VIP RCN1',
                            'Inh L1 VIP PRSS8 or TNFAIP8L3',
                            'Inh L1 ADARB2 ADAM33 or SST CXCL14',
                            'Inh L1 ADARB2 DISP2 or VIP SOX11',
                            'Inh L1-6 VIP PENK or L1-5 VIP KCNJ2',
                            'Inh L1 VIP PCDH20 or L1-2 VIP PPAPDC1A',
                            'Inh L2-6 VIP VIP or L3-6 VIP KCTD13',
                            'Inh L1-2 VIP RPL41P3 or L1-3 VIP ACHE',
                            'Inh L2-4 VIP LGI2 or L1-4 VIP CHRNA2',
                            'Inh L2-5 VIP TOX2 or L2-4 VIP DSEL',
                            'Inh L1-3 VIP ZNF322P1 or L3 VIP CBLN1',
                            'Inh L5-6 SST ISOC1 or KLHL14',
                            'Inh L3-5 SST MAFB or L4-6 SST MTHFD2P6',
                            'Inh L6 LHX6 GLP1R or L5-6 PVALB FAM150B',
                            'Inh L5 PVALB CNTNAP3P2 or L2-4 PVALB C8orf4',
                            'Inh L1-2 PVALB TAC1 or L1-3 PVALB WFDC2',
                            'Inh L1-6 PVALB SCUBE3 or L3-6 PVALB MFI2',
                            'Exc L2-4 RORB GRIK1 or L3-4 RORB RPS3P6',
                            'Exc L4 RORB BHLHE22 or CCDC168',
                            'Exc L4-5 RORB ASCL1 or AIM2',
                            'Exc L3 LINC00507 PSRC1 or L3-5 LINC00507 SLN',
                            'Exc L3 RORB CARTPT or L3-4 RORB FOLH1B',
                            'Exc L3-4 RORB SEMA6D or PRSS12',
                            'Exc L3 LINC00507 CTXN3 or THEMIS PLA2G7',
                            'Exc L5-6 THEMIS OR1J1 or L6 THEMIS EGR3',
                            'Exc L6 THEMIS LINC00343 or L5-6 THEMIS TMEM233',
                            'Exc L3-5 RORB CMAHP or CD24',
                            'Exc L5-6 RORB or L5 RORB LINC01202',
                            'Exc L4-5 RORB HNRNPA1P46 or RPL31P31',
                            'Exc L3-5 THEMIS UBE2F or L5 RORB SNHG7',
                            'Exc L6 THEMIS C6orf48 or L5-6 THEMIS GPR21',
                            'Exc L5-6 FEZF2 ANKRD20A1 or L6 FEZF2 FAM95C',
                            'Exc L6 FEZF2 TBC1D26 or KRT17',
                            'Exc L6 FEZF2 SLITRK6 or TBCC',
                            'Exc L5 FEZF2 SCN7A or MORN2',
                            'Exc L5 FEZF2 DYRK2 or L5-6 FEZF2 MYBPHL',
                            'Exc L5-6 FEZF2 CYP26B1 or RSAD2',
                            'Astro L1 FGFR3 MT1G or FOS',
                            'Oligo L4-6 MOBP COL18A1 or OPALIN',
                            'VLMC L1-3 CYP1B1 or Peri L1-6 MUSTN1']
            for k in range(0,len(new_names)):
                metadata['cell_type_alias_label_76CTs'].iloc[(metadata['cell_type_alias_label_76CTs']==old_names[k][0]) | (metadata['cell_type_alias_label_76CTs']==old_names[k][1])] = new_names[k] 
            if n_CTs == 51:
                metadata['cell_type_alias_label_51CTs'] = metadata['cell_type_alias_label_76CTs'].copy()
                old_names = [['Inh L1-4 LAMP5 DUSP4 or L6 LAMP5 C1QL2','Inh L1-6 LAMP5 CA13 or L5-6 LAMP5 SFTA3'],
                                ['Inh L1 PAX6 CA4 or GRIP2','Inh L1-3 PAX6 NABP1 or L1-6 VIP RCN1'],
                                ['Inh L1 ADARB2 ADAM33 or SST CXCL14','Inh L1 ADARB2 DISP2 or VIP SOX11'],
                                ['Inh L2-6 VIP VIP or L3-6 VIP KCTD13','Inh L1-6 VIP RGS16'],
                                ['Inh L1-3 VIP SSTR1','Inh L1-2 VIP RPL41P3 or L1-3 VIP ACHE'],
                                ['Inh L2-4 VIP LGI2 or L1-4 VIP CHRNA2','Inh L1-3 VIP CCDC184'],
                                ['Inh L1-3 VIP ZNF322P1 or L3 VIP CBLN1','Inh L1-3 VIP GGH'],
                                ['Inh L5-6 SST ISOC1 or KLHL14','Inh L3-5 SST MAFB or L4-6 SST MTHFD2P6'],
                                ['Inh L5 PVALB CNTNAP3P2 or L2-4 PVALB C8orf4','Inh L5-6 PVALB STON2'],
                                ['Exc L2-4 RORB GRIK1 or L3-4 RORB RPS3P6','Exc L4 RORB BHLHE22 or CCDC168'],
                                ['Exc L4-5 RORB ASCL1 or AIM2','Exc L4 RORB CACNG5'],
                                ['Exc L3 LINC00507 PSRC1 or L3-5 LINC00507 SLN','Exc L2-3 LINC00507 RPL9P17'],
                                ['Exc L3-4 RORB SEMA6D or PRSS12','Exc L3-5 RORB HSPB3'],
                                ['Exc L3-5 THEMIS ELOF1','Exc L3 LINC00507 CTXN3 or THEMIS PLA2G7'],
                                ['Exc L5-6 THEMIS OR1J1 or L6 THEMIS EGR3','Exc L6 THEMIS LINC00343 or L5-6 THEMIS TMEM233'],
                                ['Exc L4-6 RORB HPCA','Exc L5-6 RORB or L5 RORB LINC01202'],
                                ['Exc L3-5 THEMIS UBE2F or L5 RORB SNHG7','Exc L4-5 RORB LCN15'],
                                ['Exc L6 THEMIS C6orf48 or L5-6 THEMIS GPR21','Exc L5-6 THEMIS THTPA'],
                                ['Exc L5-6 FEZF2 ANKRD20A1 or L6 FEZF2 FAM95C','Exc L6 FEZF2 CPZ'],
                                ['Exc L6 FEZF2 ETV4','Exc L6 FEZF2 TBC1D26 or KRT17'],
                                ['Exc L6 FEZF2 P4HA3','Exc L6 FEZF2 SLITRK6 or TBCC'],
                                ['Exc L3-5 FEZF2 DCN','Exc L5 FEZF2 SCN7A or MORN2'],
                                ['Exc L5 FEZF2 DYRK2 or L5-6 FEZF2 MYBPHL','Exc L5-6 FEZF2 CYP26B1 or RSAD2'],
                                ['Astro L1-6 FGFR3 ETNPPL','Astro L1 FGFR3 MT1G or FOS'],
                                ['Endo L2-5 CLDN5','VLMC L1-3 CYP1B1 or Peri L1-6 MUSTN1']]
                new_names = ['Inh L1-4 LAMP5 DUSP4 or L6 LAMP5 C1QL2 or L1-6 LAMP5 CA13 or L5-6 LAMP5 SFTA3',
                                'Inh L1 PAX6 CA4 or GRIP2 or L1-3 PAX6 NABP1 or L1-6 VIP RCN1',
                                'Inh L1 ADARB2 ADAM33 or SST CXCL14 or ADARB2 DISP2 or VIP SOX11',
                                'Inh L2-6 VIP VIP or L3-6 VIP KCTD13 or L1-6 VIP RGS16',
                                'Inh L1-3 VIP SSTR1 or ACHE or L1-2 VIP RPL41P3',
                                'Inh L2-4 VIP LGI2 or L1-4 VIP CHRNA2 or L1-3 VIP CCDC184',
                                'Inh L1-3 VIP ZNF322P1 or GGH or L3 VIP CBLN1',
                                'Inh L5-6 SST ISOC1 or KLHL14 or L3-5 SST MAFB or L4-6 SST MTHFD2P6',
                                'Inh L5 PVALB CNTNAP3P2 or L2-4 PVALB C8orf4 or L5-6 PVALB STON2',
                                'Exc L2-4 RORB GRIK1 or L3-4 RORB RPS3P6 or L4 RORB BHLHE22 or CCDC168',
                                'Exc L4-5 RORB ASCL1 or AIM2 or L4 RORB CACNG5',
                                'Exc L3 LINC00507 PSRC1 or L3-5 LINC00507 SLN or L2-3 LINC00507 RPL9P17',
                                'Exc L3-4 RORB SEMA6D or PRSS12 or L3-5 RORB HSPB3',
                                'Exc L3-5 THEMIS ELOF1 or L3 LINC00507 CTXN3 or THEMIS PLA2G7',
                                'Exc L5-6 THEMIS OR1J1 or TMEM233 or L6 THEMIS EGR3 or LINC00343',
                                'Exc L4-6 RORB HPCA or L5-6 RORB or L5 RORB LINC01202',
                                'Exc L3-5 THEMIS UBE2F or L5 RORB SNHG7 or L4-5 RORB LCN15',
                                'Exc L6 THEMIS C6orf48 or L5-6 THEMIS GPR21 or THTPA',
                                'Exc L5-6 FEZF2 ANKRD20A1 or L6 FEZF2 FAM95C or FEZF2 CPZ',
                                'Exc L6 FEZF2 ETV4 or TBC1D26 or KRT17',
                                'Exc L6 FEZF2 P4HA3 or SLITRK6 or TBCC',
                                'Exc L3-5 FEZF2 DCN or L5 FEZF2 SCN7A or MORN2',
                                'Exc L5 FEZF2 DYRK2 or L5-6 FEZF2 MYBPHL or CYP26B1 or RSAD2',
                                'Astro L1-6 FGFR3 ETNPPL or L1 FGFR3 MT1G or FOS',
                                'Endo L2-5 CLDN5 or VLMC L1-3 CYP1B1 or Peri L1-6 MUSTN1']
                for k in range(0,len(new_names)):
                    metadata['cell_type_alias_label_51CTs'].iloc[(metadata['cell_type_alias_label_51CTs']==old_names[k][0]) | (metadata['cell_type_alias_label_51CTs']==old_names[k][1])] = new_names[k] 
        else:
            metadata['cell_type_alias_label_18CTs'] = metadata['cell_type_alias_label_18CTs'].copy()
            old_names = [['Inh L1 LAMP5 GGT8P','Inh L1 LAMP5 NDNF','Inh L1-4 LAMP5 DUSP4','Inh L6 LAMP5 C1QL2','Inh L1-6 LAMP5 CA13', 'Inh L5-6 LAMP5 SFTA3','Inh L6 LAMP5 ANKRD20A11P'],
                            ['Inh L1 PAX6 CA4','Inh L1 PAX6 GRIP2','Inh L1-3 PAX6 NABP1','Inh L1-2 PAX6 SCGN'],
                            ['Inh L1-6 VIP RCN1','Inh L1 VIP PRSS8','Inh L1 VIP TNFAIP8L3','Inh L1 VIP SOX11','Inh L1-6 VIP PENK','Inh L1-5 VIP KCNJ2','Inh L1 VIP PCDH20','Inh L1-2 VIP PPAPDC1A','Inh L2-6 VIP VIP','Inh L3-6 VIP KCTD13','Inh L1-2 VIP RPL41P3','Inh L1-3 VIP ACHE','Inh L2-4 VIP LGI2','Inh L1-4 VIP CHRNA2','Inh L2-5 VIP TOX2','Inh L2-4 VIP DSEL','Inh L1-3 VIP ZNF322P1','Inh L3 VIP CBLN1','Inh L1-6 VIP RGS16','Inh L1-3 VIP SSTR1', 'Inh L1-3 VIP CCDC184','Inh L1-3 VIP GGH'],
                            ['Inh L1 ADARB2 ADAM33','Inh L1 ADARB2 DISP2'],
                            ['Inh L1 SST CXCL14','Inh L5-6 SST ISOC1','Inh L5-6 SST KLHL14','Inh L3-5 SST MAFB','Inh L4-6 SST MTHFD2P6','Inh L6 SST NPY','Inh L5-6 SST TH','Inh L2-4 SST AHR'],
                            ['Inh L6 LHX6 GLP1R'],
                            ['Inh L5-6 PVALB FAM150B','Inh L5 PVALB CNTNAP3P2','Inh L2-4 PVALB C8orf4','Inh L1-2 PVALB TAC1','Inh L1-3 PVALB WFDC2','Inh L1-6 PVALB SCUBE3','Inh L3-6 PVALB MFI2','Inh L4-5 PVALB TRIM67','Inh L5-6 PVALB STON2','Inh L3-4 PVALB HOMER3'],
                            ['Exc L2-4 RORB GRIK1','Exc L3-4 RORB RPS3P6','Exc L4 RORB BHLHE22','Exc L4 RORB CCDC168','Exc L4-5 RORB ASCL1','Exc L4-5 RORB AIM2','Exc L3 RORB CARTPT','Exc L3-4 RORB FOLH1B','Exc L3-4 RORB SEMA6D','Exc L3-4 RORB PRSS12','Exc L3-5 RORB CMAHP','Exc L3-5 RORB CD24','Exc L5-6 RORB LINC00320','Exc L5 RORB LINC01202','Exc L4-5 RORB HNRNPA1P46','Exc L4-5 RORB RPL31P31','Exc L5 RORB SNHG7','Exc L4 RORB CACNG5','Exc L3-5 RORB HSPB3','Exc L4-6 RORB HPCA','Exc L4-5 RORB LCN15'],
                            ['Exc L3 LINC00507 PSRC1','Exc L3-5 LINC00507 SLN','Exc L3 LINC00507 CTXN3','Exc L2-3 LINC00507 RPL9P17','Exc L4-5 RORB LINC01474'],
                            ['Exc L3 THEMIS PLA2G7','Exc L5-6 THEMIS OR1J1','Exc L6 THEMIS EGR3','Exc L6 THEMIS LINC00343','Exc L5-6 THEMIS TMEM233','Exc L3-5 THEMIS UBE2F','Exc L6 THEMIS C6orf48','Exc L5-6 THEMIS GPR21','Exc L3-5 THEMIS ELOF1','Exc L5-6 THEMIS IL7R','Exc L5-6 THEMIS THTPA'],
                            ['Exc L5-6 FEZF2 ANKRD20A1','Exc L6 FEZF2 FAM95C','Exc L6 FEZF2 TBC1D26','Exc L6 FEZF2 KRT17','Exc L6 FEZF2 SLITRK6','Exc L6 FEZF2 TBCC','Exc L5 FEZF2 SCN7A','Exc L5 FEZF2 MORN2','Exc L5 FEZF2 DYRK2','Exc L5-6 FEZF2 MYBPHL','Exc L5-6 FEZF2 CYP26B1','Exc L5-6 FEZF2 RSAD2','Exc L6 FEZF2 VWA2','Exc L6 FEZF2 CPZ','Exc L6 FEZF2 ETV4','Exc L6 FEZF2 P4HA3','Exc L3-5 FEZF2 ONECUT1','Exc L3-5 FEZF2 DCN','Exc L5-6 FEZF2 CABP7'],
                            ['Astro L1 FGFR3 MT1G','Astro L1 FGFR3 FOS','Astro L1-6 FGFR3 ETNPPL'],
                            ['Oligo L4-6 MOBP COL18A1','Oligo L4-6 OPALIN']]
            new_names = ['Inh LAMP5',
                            'Inh PAX6',
                            'Inh VIP',
                            'Inh L1 ADARB2',
                            'Inh SST',
                            'Inh L6 LHX6',
                            'Inh PVALB',
                            'Exc RORB',
                            'Exc LINC00507',
                            'Exc THEMIS',
                            'Exc FEZF2',
                            'Astro FGFR3',
                            'Oligo L4-6']
            for k in range(0,len(new_names)):
                for l in range(0,len(old_names[k])):
                    metadata['cell_type_alias_label_18CTs'].iloc[(metadata['cell_type_alias_label_18CTs']==old_names[k][l])] = new_names[k] 
    else:
        print('The specified ref data set is not available!')

    metadata.rename(columns = {'sample_name':'sample', 'cell_type_alias_label_'+str(n_CTs)+'CTs':'cell type'},inplace=True)
    return metadata
    
def integrate_conos_cluster_in_metadata(df_c,path_loom_file,loom_file_name, opt_merge_samples,opt_pagoda_filtering,opt_conos_resolution,samples_to_skip):
    add_str = '_cellranger'

    if opt_pagoda_filtering:
        add_str_pagoda = '_pagoda'
    else:
        add_str_pagoda = ''
    os.chdir(path_loom_file)
    if os.path.isfile(loom_file_name):
        with loompy.connect(loom_file_name) as D:
            sample_IDs,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
    else: 
        with loompy.connect('Samples'+add_str+add_str_pagoda+'_TH_and_D_adj_filtered_and_CT_annotated.loom') as D:
            sample_IDs,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
    counter=0
    for s_id,samplenr in enumerate(sample_IDs):
        if samplenr in samples_to_skip:
            continue
        if opt_pagoda_filtering:
            loom_file_name_i = samplenr + loom_file_name[7:-22]+'.loom'
        else:
            loom_file_name_i = samplenr + loom_file_name[7:]
        if s_id==0:
            target_files = [loom_file_name_i]
        else:
            target_files.append(loom_file_name_i)
        Di = loompy.connect(target_files[counter])
        CellIDs_interest = Di.ca.CellID
        Conos_Di = [df_c[df_c['Cell_ID']==cID]['conos_cluster'].values[0] if df_c[df_c['Cell_ID']==cID]['conos_cluster'].values.size !=0 else 'removed' for cID in CellIDs_interest]
        Conos_Di_A = np.array(Conos_Di)
        #keep existing row and column attributes as they are but add/ replace the one to be integrated
        add_col_metadata_to_loom_file(Di, Conos_Di_A, "Conos_cluster_res_"+opt_conos_resolution, path_loom_file, target_files[counter], False)
        print(samplenr+' done!')
        counter=counter+1
    if opt_merge_samples==True:
        os.chdir(path_loom_file)
        loompy.combine(files = target_files, output_file = loom_file_name, batch_size = 3000)

def add_col_metadata_to_loom_file(D, new_col_attr, new_col_attr_name, file_path, filename, opt_automatically_adapt_file_name):
    row_attrs, col_attrs = get_col_and_row_attributes_as_dictionary(D)
    #col_attrs = D.ca
    col_attrs[new_col_attr_name]=new_col_attr
    os.chdir(file_path)
    M = D[:,:]
    if opt_automatically_adapt_file_name:
        new_filename = filename[:-5]+'_with_'+new_col_attr_name+'.loom'
    else:
        new_filename = filename
    D.close()
    loompy.create(new_filename, M, row_attrs, col_attrs)

    return new_filename

def get_col_and_row_attributes_as_dictionary(D):
    bool_first = True
    for rk in D.ra.keys():
        if bool_first:
            row_attrs = {rk: D.ra[rk]}
            bool_first = False
        else:
            row_attrs[rk] = D.ra[rk]
    #row_attrs = D.ra
    bool_first = True
    for ck in D.ca.keys():
        if bool_first:
            col_attrs = {ck: D.ca[ck]}
            bool_first = False
        else:
            col_attrs[ck] = D.ca[ck]
    return row_attrs, col_attrs
    
def create_loom_file(D, genes_kept, cells_kept, path, filename):
    #row_attrs, col_attrs = get_col_and_row_attributes_as_dictionary_for_selection_of_cells_and_genes(D, genes_kept, cells_kept)
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
    
def create_subsample_file_percentage_cells(loom_file_path, loom_file_name, target_percentage_of_cells, opt_save_as_anndata):
    new_loom_file_name = loom_file_name[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom'
    os.chdir(loom_file_path)  
    
    if os.path.isfile(loom_file_path+"S1"+loom_file_name[7:]):
        #if file per donor exists: use those
        with loompy.connect(loom_file_name,validate=False) as D:
            #all genes will be kept
            genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
            # 20% of cells (for each sample) will be kept
            #get IDs for cell belonging to a sample --> pick target_percentage_of_cells % of IDs and add to vector cells_kept
            sample_IDs,n_cells,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
            #n_cells_target = np.round(n_cells*(target_percentage_of_cells/100)).astype(int)
        
        gc.collect()
        
        for d_id in range(0,len(sample_IDs)):
            cells_kept = []
            os.chdir(loom_file_path)  
            with loompy.connect(sample_IDs[d_id]+loom_file_name[7:]) as D:
                #all genes will be kept
                genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
                upper_val = np.shape(D)[1]
                n_cells_target = np.round(upper_val*(target_percentage_of_cells/100)).astype(int)
                cells_kept_IDs = random.sample(range(0,upper_val),n_cells_target)
                cells_kept = np.zeros((1,upper_val), dtype=bool).flatten()
                cells_kept[cells_kept_IDs] = True
                if d_id==0:
                    target_files = [sample_IDs[d_id]+loom_file_name[7:]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom']
                else:
                    target_files = target_files + [sample_IDs[d_id]+loom_file_name[7:]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom']
                #create 1 file per sample
                create_loom_file(D, genes_kept, cells_kept,loom_file_path,target_files[d_id])
            #plt.scatter(range(0,np.sum(n_cells)),range(0,np.sum(n_cells)),s=1,color='blue',marker='o') 
            #plt.scatter(cells_kept,cells_kept,s=1,color='red',marker='.')
        #merge sample specific files:
        loompy.combine(files = target_files, output_file = new_loom_file_name, batch_size = 3000)
    else: 
        #otherwise take file with all samples
        with loompy.connect(loom_file_name, validate=False) as D:
            #all genes will be kept
            genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
            # 20% of cells (for each sample) will be kept
            #get IDs for cell belonging to a sample --> pick target_percentage_of_cells % of IDs and add to vector cells_kept
            sample_IDs,n_cells,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
            #n_cells_target = np.round(n_cells*(target_percentage_of_cells/100)).astype(int)
            for d_id in range(0,len(sample_IDs)):
                #all genes will be kept
                genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
                donor_cells_idx = np.asarray(np.where(D.ca.Donor==sample_IDs[d_id]))
                upper_val_i = np.shape(donor_cells_idx)[1] #np.shape(D[:,donor_cells_idx][1]
                n_cells_target_i = np.round(upper_val_i*(target_percentage_of_cells/100)).astype(int)
                cells_kept_IDs_i = donor_cells_idx[:,random.sample(range(0,upper_val_i),n_cells_target_i)]#[donor_cells_idx[i] for i in random.sample(range(0,upper_val_i),n_cells_target_i)]
                if d_id==0:
                    cells_kept_IDs = cells_kept_IDs_i[0][:]
                else:
                    cells_kept_IDs = np.append(cells_kept_IDs,cells_kept_IDs_i[0][:])
            cells_kept = np.zeros((1,np.shape(D)[1]), dtype=bool).flatten()
            cells_kept[cells_kept_IDs] = True
            create_loom_file(D, genes_kept, cells_kept,loom_file_path,new_loom_file_name)
    if opt_save_as_anndata:
        aD = anndata.read_loom(loom_file_name[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom', X_name='matrix') 
        aD.var_names_make_unique()
        aD.write_h5ad(filename = loom_file_name[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.h5ad')
        
    return new_loom_file_name

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
    df = pd.read_excel(path_filtered_data+'Cell_type_colors.xlsx',engine='openpyxl')
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

def plot_Sankey_of_CT_clusters(CT_class,D,ca_variable,n_clusters,opt_sankey_without_removed_cells,opt_sankey_without_unassigned_cells,opt_sankey_of_merged_cluster,df_ct_mapping,opt_D_version,path_results_figures):
    if opt_sankey_of_merged_cluster:
        DF = pd.DataFrame({'Cluster':D.ca['cluster_'+str(n_clusters)+'CTs'], 'Celltype':D.ca[ca_variable+'_red']})
    else:
        DF = pd.DataFrame({'Cluster':D.ca['Conos_2nd_level'], 'Celltype':D.ca[ca_variable]})
    if opt_sankey_without_unassigned_cells:
        DF = DF[DF.Cluster!='unassigned']
    if opt_sankey_without_removed_cells:
        DF = DF[DF.Cluster!='removed']
        DF = DF[DF.Cluster!='nan']
        DF = DF[DF.Cluster!='other']
    F = pd.crosstab(DF['Cluster'],DF['Celltype']).sort_values('Cluster', ascending=True)
    if opt_sankey_of_merged_cluster:
        bool_first=True
        for s in F.index.tolist():
            if bool_first:
                colors = df_ct_mapping[df_ct_mapping["cluster"]==s]["color"].tolist()
                bool_first=False
            else:
                colors = colors + df_ct_mapping[df_ct_mapping["cluster"]==s]["color"].tolist()
        F['colors_source'] = colors
    else:    
        bool_first=True
        for s in F.index.tolist():
            if bool_first:
                colors = df_ct_mapping[df_ct_mapping["cluster"]==s]["color"].tolist()
                bool_first=False
            else:
                colors = colors + df_ct_mapping[df_ct_mapping["cluster"]==s]["color"].tolist()
        F['colors_source'] = colors
        #F['colors_source'] = ['#CC0000']*len(F.index.tolist())

    n_groups_left = len(F.index.tolist())
    n_groups_right = len(F.columns.tolist())-1
    label = np.concatenate([F.index.tolist(),F.columns[:-1].tolist()])
    source = []
    for id_left in range(0,n_groups_left):
        source = source + [id_left]*n_groups_right #source nodes (left part of the plot)
    ids_right = np.arange(n_groups_left,n_groups_left+n_groups_right)
    target = np.tile(ids_right,n_groups_left).tolist()
    value=[]
    for cl in F.index.tolist():
        value = value + F[F.index==cl].values[0][:-1].tolist()
    link = dict(source = source, target = target, value=value)#, color=color_link)
    node = dict(label = label, pad=55, line = dict(color = "black", width = 0.5), thickness=15, color = F['colors_source']) # pad = space between categories in the left
    data = go.Sankey(link = link, node=node)
    fig = go.Figure(data)
    fig.update_layout(
        #hovermode = 'x',
        title = CT_class+" Cells",
        #paper_bgcolor = '#51504f',
        font = dict(size=15,color='black')  
    )
    fig.show()
    #opens it in a browser:
    #plot(fig, auto_open=True)
    if opt_sankey_without_removed_cells:
        add_str_sankey = '_without_removed_cells'
    else:
        add_str_sankey = ''
    fig.write_html(path_results_figures+"figures_"+opt_D_version+"_Sankey/"+CT_class+'_cluster_'+str(n_clusters)+'CTs'+ca_variable[14:-1]+add_str_sankey+".html")
    #fig.show()

def plot_PCA_UMAP_clusterContribution_from_loom(loom_file_name, loom_file_path, add_str_figure_file_name,path_results,variables_color_code,single_marker_gene_list,single_marker_CT_list,lei_res,opt_server,opt_clustering_plots_single_sample,opt_which_clustering,n_clusters_range,opt_abundance,opt_sep_4_scz_ctrl,opt_plot_umap):
    if opt_server:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt 
        
    #get aD in right shape:
        #-attributes as type str
        #-remove xy chromosome genes
        #-add colors for col attributes
        #-normalization
        #-scaling
        #-HVGs
        #-PCA
        #-neighborhood graph for UMAP
        #-UMAP calculation
    aD = get_anndata_ready4plotting(loom_file_path,loom_file_name,opt_which_clustering,n_clusters_range,variables_color_code)
    
    
    os.chdir(path_results)
    if len(single_marker_gene_list)>0:
        if 'CT_class_conos' in variables_color_code:
            sc.pl.stacked_violin(aD, single_marker_gene_list, groupby='CT_class_conos',swap_axes=True,save="_violin_CT_class_conos.pdf")
        if 'CT_class_scmap' in variables_color_code:
            sc.pl.stacked_violin(aD, single_marker_gene_list, groupby='CT_class_scmap',swap_axes=True,save="_violin_CT_class_scmap.pdf")
        
    if opt_clustering_plots_single_sample:
        metadata_variables=['CT_ann_ABM_MCA_scmap_cluster']
    else:
        metadata_variables = ['Disease', 'Donor', 'Library', 'Sex']
        if ('cluster_name_15CTs' in variables_color_code):
            metadata_variables = metadata_variables+['cluster_name_15CTs']
        if ('CT_ann_ABM_MCA_scmap_cluster_51_CTs' in variables_color_code):
            metadata_variables = metadata_variables + ['CT_ann_ABM_MCA_scmap_cluster_51_CTs']
        if ('CT_ann_ABM_MCA_scmap_cluster_76_CTs' in variables_color_code):
            metadata_variables=metadata_variables +['CT_ann_ABM_MCA_scmap_cluster_76_CTs']
        if ('CT_ann_ABM_MCA_scmap_cell2cluster_51CTs' in variables_color_code):
            metadata_variables=metadata_variables +['CT_ann_ABM_MCA_scmap_cell2cluster_51_CTs']
        if ('CT_ann_ABM_MCA_scmap_cell2cluster_76CTs' in variables_color_code):
            metadata_variables=metadata_variables +['CT_ann_ABM_MCA_scmap_cell2cluster_76_CTs']

    n_cells = np.shape(aD)[0]                   
    for v_id, v in enumerate(variables_color_code):
        if v.startswith('cluster_name_') & v.endswith('CTs'):
            n_clusters = int(v.split('cluster_name_')[1].split('CTs')[0])
            
            #remove the removed cells
            bool_cells_kept = aD.obs['cluster_name_'+str(n_clusters)+'CTs']!='none (removed)'
            aD_red = aD[bool_cells_kept,:]
            
            CT_order_list = get_unique_values_in_list_while_preserving_order(aD.obs['cluster_name_'+str(n_clusters)+'CTs'].tolist())
            CT_order_list_red = get_unique_values_in_list_while_preserving_order(aD_red.obs['cluster_name_'+str(n_clusters)+'CTs'].tolist())
            
            #reorder cells according to cell types:
            # add colors for 15 and 37 CTs to aD_red
            CT_order_list_red_reordered = np.array(CT_order_list_red)[np.char.startswith(CT_order_list_red,'Exc')].tolist() + np.array(CT_order_list_red)[np.char.startswith(CT_order_list_red,'Inh')].tolist() + np.array(CT_order_list_red)[np.logical_and(~np.char.startswith(CT_order_list_red,'Exc'),~np.char.startswith(CT_order_list_red,'Inh'))].tolist()
            aD_red.obs['cluster_name_'+str(n_clusters)+'CTs'] = aD_red.obs['cluster_name_'+str(n_clusters)+'CTs'].astype('category')
            aD_red.obs['cluster_name_'+str(n_clusters)+'CTs'] = aD_red.obs['cluster_name_'+str(n_clusters)+'CTs'].cat.reorder_categories(CT_order_list_red_reordered)
            aD_red.uns['cluster_name_'+str(n_clusters)+'CTs_colors'] = get_colors_for_CTs(CT_order_list_red_reordered,n_clusters,aD_red.obs)
            
            n_cells_red = np.shape(aD_red)[0]
            
            #plot umaps:
            if opt_plot_umap:
                #plot UMAP
                sc.pl.umap(aD_red, color=v, size=100000/n_cells_red, save="_"+add_str_figure_file_name+v+".pdf")#,
                print('done deal UMAP!')
                #sc.pl.umap(aD, color=v, color_map = np.unique(aD.obs["cluster_15CTs_colors"]).tolist(), palette = np.unique(aD.obs["cluster_15CTs_colors"]).tolist(), save="_"+v+".pdf")#,
            
            #plot violin of number of genes with and without removed cells:
            plot_violins_for_number_of_genes(aD, n_clusters, CT_order_list, False, add_str_figure_file_name+"_with_removed")
            plot_violins_for_number_of_genes(aD_red, n_clusters, CT_order_list_red, True, add_str_figure_file_name)    
        
            #plot dotplot marker genes expression per CT:
            sc.tl.dendrogram(aD_red,groupby='cluster_name_'+str(n_clusters)+'CTs')
            marker_genes_dict = {'Excitatory neurons ': ['SLC17A7','CBLN4','CUX2','RORB','THEMIS','FEZF2'],
                                'Inhibitory neurons': ['GAD1','TAC1','PVALB','VIP','CHRNA2','CALB2','LAMP5','SV2C','SST', 'LHX6'],
                                'Astrocytes':['AQP4'],
                                'Microglial cells':['CSF1R'],
                                'Oligodendrocyes':['MOG'],
                                'Oligodendrocyte progenitor cells':['PDGFRA']
                                }
            
            marker_genes_Ruzicka_dict = {'Excitatory neurons ': ['SLC17A7', 'NRGN'], 
                                         'Excitatory L2/3 neurons ': ['CUX2'],
                                         'Excitatory L4 neurons ': ['RORB'],
                                         'Excitatory L5/6 neurons ': ['TELE4'], 
                                         'Inhibitory neurons': ['GAD2'],
                                         'Inhibitory Reelin': ["RELN"],
                                         'Inhibitory PVALB': ["PV", "PVALB"],
                                         'Inhibitory SST': ["SST"],
                                         'Inhibitory Rosehip': ["SV2C"],
                                         'Astrocytes':['SLC1A2'],
                                         'Oligodendrocyes':['PLP1'],
                                         'Oligodendrocyte progenitor cells':['VCAN'],
                                         'Microglial cells':['CSF1R']}
            
            sc.pl.dotplot(aD_red, marker_genes_dict, 'cluster_name_'+str(n_clusters)+'CTs', dendrogram=False, categories_order=CT_order_list_red, save=add_str_figure_file_name+'_cluster_name_'+str(n_clusters)+'CTs')
            sc.pl.dotplot(aD_red, marker_genes_Ruzicka_dict, 'cluster_name_'+str(n_clusters)+'CTs', dendrogram=False, categories_order=CT_order_list_red, save=add_str_figure_file_name+'_cluster_name_'+str(n_clusters)+'CTs_Ruzicka')
            print('done deal violin plots!')

            plot_abundance_of_cells_from_each_condition_per_cluster(aD,metadata_variables,v,opt_server,path_results,add_str_figure_file_name,opt_abundance,opt_sep_4_scz_ctrl)
        elif (v=='CT_ann_ABM_MCA_scmap_cluster_51CTs') or (v =='CT_ann_ABM_MCA_scmap_cluster_76CTs') and opt_plot_umap:
            sc.pl.umap(aD, color=v, size=100000/n_cells, color_map = aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"], palette = aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"], save="_"+add_str_figure_file_name+v+".pdf")#, palette="twilight", color_map=plt.cm.Reds)
        else:
            if opt_plot_umap:
                sc.pl.umap(aD, color=v, size=100000/n_cells, save="_"+add_str_figure_file_name+v+".pdf")#, palette="twilight", color_map=plt.cm.Reds)
    

    for r in lei_res:
        clustering_key="leiden_"+str(r)
        sc.tl.leiden(aD,resolution = r,key_added = clustering_key)
        if opt_plot_umap:
            sc.pl.umap(aD, color=[clustering_key],save="_"+add_str_figure_file_name+clustering_key, size=100000/n_cells)
            plt.close()
        
        #plot UMAP and heatmap
        #Plot proportion of cells from each condition (Disease, Sex or annotated Cell type, Donor, Library) per cluster:
        plot_abundance_of_cells_from_each_condition_per_cluster(aD,metadata_variables,clustering_key,opt_server,path_results,add_str_figure_file_name,opt_abundance,opt_sep_4_scz_ctrl)

    if len(np.unique(aD.obs['Donor']))>1:     
        if opt_plot_umap:
            for m_id,m in enumerate(single_marker_gene_list):
                sc.pl.umap(aD, color=m, size=100000/n_cells, save="_"+add_str_figure_file_name+single_marker_CT_list[m_id]+".pdf")#, palette="twilight", color_map=plt.cm.Reds)
    #os.chdir(loom_file_path)
    #del aD.layers['spliced']
    #del aD.layers['unspliced']
    #aD.write(filename=add_str_figure_file_name+".h5ad")

def get_unique_values_in_list_while_preserving_order(input_list):
    _, idx = np.unique(input_list, return_index=True)
    idx_sorted = sorted(idx)
    unique_list_with_preserved_order = [input_list[i] for i in idx_sorted]
    
    return unique_list_with_preserved_order


def plot_violins_for_number_of_genes(aD, n_clusters, CT_order_list, opt_red, add_str_figure_file_name):
    # mitochondrial genes
    aD.var['mt'] = aD.var_names.str.startswith('MT-') 
    sc.pp.calculate_qc_metrics(aD, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    aD.obs['pct_counts_mt'] = np.sum(aD[:,aD.var['mt']==True].X.toarray()!=0,1)*100/np.sum(aD.X.toarray()!=0,1)
    # add the total counts per cell as observations-annotation to adata
    aD.obs['n_genes'] = np.sum(aD.X.toarray() != 0, 1)
    # CT_names_reordered_w_r = CT_names_reordered
    # if opt_red == False:
    #     CT_names_reordered_w_r = CT_names_reordered_w_r.append('none (removed)')
    sc.pl.violin(aD, ['n_genes', 'pct_counts_mt'], groupby='cluster_name_'+str(n_clusters)+'CTs', categories_order=CT_order_list, stripplot=False, rotation=90,  inner='box',save="_"+add_str_figure_file_name+'cluster_name_'+str(n_clusters)+"CTs.pdf")  # use stripplot=False to remove the internal dots, inner='box' adds a boxplot inside violins


def plot_abundance_of_cells_from_each_condition_per_cluster(aD,metadata_variables, clustering_key,opt_server,path_figure,file_name_str,opt_abundance,opt_sep_4_scz_ctrl):
    if opt_abundance=='relative':
        fig_str = 'P_'
    else:
        fig_str = 'N_'
    if opt_server:
        import matplotlib
        matplotlib.use('pdf')
    import matplotlib.pyplot as plt 
    pos = [1, -0.15]
    aD = get_color_maps_for_anndata(aD)
    if opt_sep_4_scz_ctrl==True:
        aD_SCZ = aD[aD.obs["Disease"]=="SCZ",:]
        aD_CTRL = aD[aD.obs["Disease"]=="CTRL",:]
    for metadata_var in metadata_variables:
        if opt_sep_4_scz_ctrl==True and metadata_var!="Disease":
            if opt_abundance=="relative":
                tmp_SCZ = pd.crosstab(aD_SCZ.obs[clustering_key],aD_SCZ.obs[metadata_var],normalize='index')
                tmp_CTRL = pd.crosstab(aD_CTRL.obs[clustering_key],aD_CTRL.obs[metadata_var],normalize='index')
            else:
                tmp_SCZ = pd.crosstab(aD_SCZ.obs[clustering_key],aD_SCZ.obs[metadata_var])
                tmp_CTRL = pd.crosstab(aD_CTRL.obs[clustering_key],aD_CTRL.obs[metadata_var])
        else:
            if opt_abundance=="relative":
                tmp = pd.crosstab(aD.obs[clustering_key],aD.obs[metadata_var],normalize='index')
            else:
                tmp = pd.crosstab(aD.obs[clustering_key],aD.obs[metadata_var])
        if (metadata_var=='CT_ann_ABM_MCA_scmap_cluster_76CTs') or (metadata_var=='CT_ann_ABM_MCA_scmap_cluster_51CTs') or (metadata_var=='CT_ann_ABM_MCA_scmap_cell2cluster_76CTs') or (metadata_var=='CT_ann_ABM_MCA_scmap_cell2cluster_51CTs'):
            cols=tmp.columns.tolist()
            nan_unassigned_id = [id for id,ct in enumerate(cols) if ct.startswith('nan') or ct.startswith('unassigned')]
            non_neuronal_id = [id for id,ct in enumerate(cols) if not ct.startswith('Inh') and not ct.startswith('Exc') and id not in nan_unassigned_id]
            neuronal_id = [id for id,ct in enumerate(cols) if id not in non_neuronal_id and id not in nan_unassigned_id]
            new_cols = tmp.columns[non_neuronal_id].tolist()+tmp.columns[neuronal_id].tolist()+tmp.columns[nan_unassigned_id].tolist()
            tmp = tmp[new_cols]
            col_map = aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"]
        elif metadata_var in ['Disease','Sex','cluster_name_15CTs']:
            col_map = colors.LinearSegmentedColormap.from_list("", aD.uns[metadata_var+"_colors"])
            #col_map = colors.ListedColormap(aD.uns[metadata_var+"_colors"])
        elif metadata_var=='Donor':
            col_map = 'tab20b'
        else:
            col_map = 'viridis'
        if opt_sep_4_scz_ctrl==True and metadata_var!="Disease":
            #to do:
            fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)
            #subplots with shared x axis
            #divide data in scz and ctrl
            #plot
            #ax1.plot(kind="bar",stacked=True,data=tmp_SCZ,colormap=col_map)
            tmp_SCZ.plot.bar(stacked=True,colormap=col_map,ax=ax1)
            ax1.legend(bbox_to_anchor=pos, loc='lower left')
            tmp_CTRL.plot.bar(stacked=True,colormap=col_map,ax=ax2)
            ax2.legend(bbox_to_anchor=pos, loc='lower left')
            ax2.set_xticklabels(ax2.get_xticklabels() , rotation=90, fontsize=4)
            #ax[0].set_xticklabels([])
        else:
            tmp.plot.bar(stacked=True,colormap=col_map).legend(bbox_to_anchor=pos, loc='lower left')
            plt.xticks(fontsize=4, rotation=90)
        #save:
        plt.savefig(path_figure+fig_str + 'cells_'+file_name_str+'_'+metadata_var+'_'+clustering_key+'.pdf', bbox_inches='tight')
        #plt.show()
        if (metadata_var=='CT_ann_ABM_MCA_scmap_cluster_76CTs') or (metadata_var=='CT_ann_ABM_MCA_scmap_cluster_51CTs'):
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter('P_cells_'+file_name_str+'_'+metadata_var+'_'+clustering_key+'_.xlsx', engine='xlsxwriter')
            tmp.to_excel(writer,sheet_name='raw',index=True)
            tmp.max(axis=1).to_excel(writer,sheet_name='max_per_cluster',index=True)
            tmp.max(axis=0).to_excel(writer,sheet_name='max_per_celltype',index=True)
            writer.save() 

def get_anndata_ready4plotting(loom_file_path,loom_file_name,opt_which_clustering,n_clusters_range,variables_color_code):
    os.chdir(loom_file_path)
    print(loom_file_name)
    aD_i = anndata.read_loom(loom_file_name, X_name='matrix') 
    aD_i.var_names_make_unique()
    #change type of age from int to float for plotting UMAP:
    #aD_i.obs['Age_c'] = aD_i.obs['Age'].astype('int')
    if opt_which_clustering == 'conos' or opt_which_clustering == 'conos_and_leiden' or opt_which_clustering == 'all':
        if 'Conos_cluster_res_5' in variables_color_code:
            aD_i.obs['Conos_cluster_res_5'] = aD_i.obs['Conos_cluster_res_5'].astype('str')
        if 'Conos_cluster_res_7' in variables_color_code:
            aD_i.obs['Conos_cluster_res_7'] = aD_i.obs['Conos_cluster_res_7'].astype('str')
    if opt_which_clustering == 'final_conos_based' or opt_which_clustering == 'all':
        for n_clusters in n_clusters_range:
            aD_i.obs['cluster_'+str(n_clusters)+'CTs'] = aD_i.obs['cluster_'+str(n_clusters)+'CTs'].astype('str')
    # for i in range(0,2):
    #     if i==0:
    #         add_str_chr_rm = '_with_xy_chr_info'
    #         aD = aD_i
    #     elif i==1:
    #         add_str_chr_rm = '_without_xy_chr_info'
    #         #remove x and y chromosome information:
    #         non_xy_chr_genes_list = aD_i.var_names[np.logical_and(aD_i.var['Chromosome']!='Y',aD_i.var['Chromosome']!='X')]
    #         aD = aD_i[:,non_xy_chr_genes_list]
    
    #remove x and y chromosome information:
    DF_xy_chr = pd.read_csv(loom_file_path+'XYchr_genes.csv',delimiter=',',index_col=None)
    non_xy_chr_genes_list =[g for a,g in zip(aD_i.var['Accession'],aD_i.var.index) if a not in DF_xy_chr['Accession'].tolist()]
    aD = aD_i[:,non_xy_chr_genes_list]  
    #add_str_chr_rm = '_without_xy_chr_info'
    
    #reorder cells according to cell types:
    # add colors for 15 and 37 CTs to aD
    if opt_which_clustering == 'final_conos_based' or opt_which_clustering == 'all':
        for n_clusters in n_clusters_range:
            CT_names = get_unique_values_in_list_while_preserving_order(aD.obs['cluster_name_'+str(n_clusters)+'CTs'].tolist())
            CT_names_reordered = np.array(CT_names)[np.char.startswith(CT_names,'Exc')].tolist() + np.array(CT_names)[np.char.startswith(CT_names,'Inh')].tolist() + np.array(CT_names)[np.logical_and(~np.char.startswith(CT_names,'Exc'),~np.char.startswith(CT_names,'Inh'))].tolist()
            aD.obs['cluster_name_'+str(n_clusters)+'CTs'] = aD.obs['cluster_name_'+str(n_clusters)+'CTs'].astype('category')
            aD.obs['cluster_name_'+str(n_clusters)+'CTs'] = aD.obs['cluster_name_'+str(n_clusters)+'CTs'].cat.reorder_categories(CT_names_reordered)
            aD.uns['cluster_name_'+str(n_clusters)+'CTs_colors'] = get_colors_for_CTs(CT_names_reordered,n_clusters,aD.obs)
        
    #get colors for other col attributes for plotting:
    aD = get_color_maps_for_anndata(aD) 
    
    ##basic normalization
    #normalize each cell by total counts over all genes, so that every cell has the same total count after normalization 
    #normalize to 10 000 reads per cell --> counts become comparable amongst cells

    sc.pp.normalize_total(aD, layers = 'all')
    ##log-transformation of data
    #Computes X=log(X+1), where log denotes the natural logarithm unless a different base is given.
    #chunked=True saves memory
    sc.pp.log1p(aD,chunked=True)#

    #highly variable genes:
    sc.pp.highly_variable_genes(aD,flavor='cell_ranger',n_top_genes=2000,inplace=True)
    
    #pca
    sc.tl.pca(aD, svd_solver='arpack')
    #sc.pl.pca(aD)
    #sc.pl.pca_variance_ratio(aD,save = "PCA_variance_ratio_"+add_str_chr_rm+'.pdf')
    #sc.pl.pca_loadings(aD,[1,2,3,4,5,6,7,8,9,10],save = "PCA_loadings_"+add_str_chr_rm+'.pdf')
    #variance_ratio_explained_PCA = aD.uns["pca"]["variance_ratio"]
    #X_pca = aD.obsm['X_pca']
    #np.save('PCA_var_explained_'+add_str_chr_rm,variance_ratio_explained_PCA,allow_pickle=True)
    #np.save('PCA_X_matrix_'+add_str_chr_rm,X_pca,allow_pickle=True)
    sc.pp.neighbors(aD, n_neighbors=10, n_pcs=30)
    sc.tl.umap(aD,random_state=192785)
    
    return aD

def get_color_maps_for_anndata(aD):
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_76CTs_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_51_CTs_colors"] = 'nipy_spectral'
    aD.uns["Disease_colors"] = ['#487BAF','#FA9835']#['steelblue','darkorange']
    aD.uns["Sex_colors"] = ['#CB5C93','#0080FF']#['palevioletred','dodgerblue']
    aD.uns["filtering_status_colors"] = ["#CC0066","#078383","#939393","#2284C6","#000000"]
    return aD

def get_colors_for_CTs(CT_names,n_clusters,aD_obs):
    palette_colors=[]
    for ct in CT_names:
        if ct == 'none (removed)':
            ct_color = ["#C0C0C0"] # grey
        else:
            ct_color = np.unique(aD_obs[aD_obs['cluster_name_'+str(n_clusters)+'CTs']==ct]["cluster_"+str(n_clusters)+"CTs_colors"]).tolist()
        #make sure it has the same order:
        palette_colors = palette_colors + ct_color
    #print(palette_colors)
    
    return palette_colors

def plot_percentage_unassigned_per_sample_and_method(UA_per_donor,ca_ann,sample_IDs,sample_IDs_CTRL,sample_IDs_SCZ,path_results,opt_server):
    plt = loadPltSettings(20,30)
    fig, axes = plt.subplots(len(ca_ann),2,figsize=(20,10),sharey=True,sharex=False)
    fig.subplots_adjust(hspace = .25, wspace = 0.02)  
    if len(ca_ann)>1:
        for k in range(0,len(ca_ann)):
            for g in range(0,2):
                if g==0:
                    SIDs = sample_IDs_CTRL
                    c='steelblue'
                else:
                    SIDs = sample_IDs_SCZ
                    c='darkorange'
                axes[k,g].bar(range(0,len(SIDs)),UA_per_donor[k][SIDs],color=c,edgecolor='black')
                axes[k,g].set_xticks(range(0,len(SIDs)))
                axes[k,g].set_xticklabels(sample_IDs[SIDs], rotation=60, ha='center',fontsize=14)
                if g==0:
                    axes[k,g].set_ylabel('Unassigned cells [%] \n'+ca_ann[k][7:])
                axes[k,g].set_xlabel('Sample IDs')
                axes[k,g].set_xlim(-1,len(SIDs))
    else:
        for g in range(0,2):
            if g==0:
                SIDs = sample_IDs_CTRL
                c='steelblue'
            else:
                SIDs = sample_IDs_SCZ
                c='darkorange'
            axes[g].bar(range(0,len(SIDs)),UA_per_donor[0][SIDs],color=c,edgecolor='black')
            axes[g].set_xticks(range(0,len(SIDs)))
            axes[g].set_xticklabels(sample_IDs[SIDs], rotation=60, ha='center',fontsize=14)
            if g==0:
                axes[g].set_ylabel('Unassigned cells [%] \n'+ca_ann[0][7:])
            axes[g].set_xlabel('Sample IDs')
            axes[g].set_xlim(-1,len(SIDs))
    #plt.savefig(path_results+'percentage_unassigned_per_sample_and_method.png', bbox_inches='tight')
    plt.savefig(path_results+'Perc_unassigned_per_sample_and_method.pdf', bbox_inches='tight')
    if not opt_server:
        plt.show()
    plt.close()
    plt.clf()

def plot_percentage_assigned_to_each_CT_per_method(ca_ann,DF_CT_A,path_results,opt_server):
    plt = loadPltSettings(20,30)
    plt.rcParams["legend.numpoints"] = 1
    fig, axes = plt.subplots(len(ca_ann),1,figsize=(20,20),sharey=True,sharex=True)
    fig.subplots_adjust(hspace = .1, wspace = 0.02)  
    def_colors = ['#2D71A5','#F08E2C']
    new_pal = dict(CTRL=def_colors[0],SCZ = def_colors[1])
    if len(ca_ann)>1:
        for k in range(0,len(ca_ann)):
            sns.stripplot(x="CT_assigned", y="Annotation_percentage", hue="Group", data=DF_CT_A[DF_CT_A['method']==ca_ann[k]], dodge=True, palette=new_pal, ax=axes[k])
            handles,labels = axes[k].get_legend_handles_labels()
            axes[k].legend(handles,labels,numpoints=1)
            axes[k].set_ylabel('Annotated cells [%] \n'+ca_ann[k][7:])
            if k==len(ca_ann)-1:
                axes[k].set_xlabel('Cell types')
                axes[k].set_ylim(0,np.max(DF_CT_A["Annotation_percentage"])+3)
                axes[k].set_xticklabels(axes[k].get_xticklabels(),rotation=90)
            else:
                axes[k].set_xlabel('')
    else:
        sns.stripplot(x="CT_assigned", y="Annotation_percentage", hue="Group", data=DF_CT_A[DF_CT_A['method']==ca_ann[0]], dodge=True, palette=new_pal, ax=axes)
        handles,labels = axes.get_legend_handles_labels()
        axes.legend(handles,labels,numpoints=1)
        axes.set_ylabel('Annotated cells [%] \n'+ca_ann[0][7:])
        axes.set_xlabel('Cell types')
        axes.set_ylim(0,np.max(DF_CT_A["Annotation_percentage"])+3)
        axes.set_xticklabels(axes.get_xticklabels(),rotation=90)
    #plt.savefig(path_results+'percentage_assigned_to_each_CT_per_sample_and_method.png', bbox_inches='tight')
    plt.savefig(path_results+'perc_assigned_to_each_CT_per_sample_and_method.pdf', bbox_inches='tight')
    if not opt_server:
        plt.show()
    plt.close()
    plt.clf()
    
def plot_heatmap_CT_vs_sample(DF_CT_A,ca_ann,path_results_figures,option,opt_server):
    for k in range(0,len(ca_ann)):
        loadPltSettings(4,2)
        if option=='percentage_assigned_to_each_CT_per_sample':
            data=DF_CT_A[DF_CT_A['method']==ca_ann[k]].pivot(index='Sample_ID',columns='CT_assigned', values='Annotation_percentage')
            #sns.heatmap(data,vmin=np.min(CT_A_per_donor[k,:]),vmax=np.max(CT_A_per_donor[k,:]),center=np.mean(CT_A_per_donor[k,:]))
            sns.heatmap(data,vmin=0,vmax=85,center=25)
        elif option=='mean_score_assigned_to_each_CT_per_sample':
            data=DF_CT_A[DF_CT_A['method']==ca_ann[k]].pivot(index='Sample_ID',columns='CT_assigned', values='mean_score')
            #sns.heatmap(data,vmin=np.min(CT_A_per_donor[k,:]),vmax=np.max(CT_A_per_donor[k,:]),center=np.mean(CT_A_per_donor[k,:]))
            sns.heatmap(data)
        elif option=='median_score_assigned_to_each_CT_per_sample':
            data=DF_CT_A[DF_CT_A['method']==ca_ann[k]].pivot(index='Sample_ID',columns='CT_assigned', values='median_score')
            #sns.heatmap(data,vmin=np.min(CT_A_per_donor[k,:]),vmax=np.max(CT_A_per_donor[k,:]),center=np.mean(CT_A_per_donor[k,:]))
            sns.heatmap(data)
        else:
            print('This option is not implemented. Go to self_fun_plot.py and add option in function plot_heatmap_CT_vs_sample!')  
        #plt.savefig(path_results_figures+option+'_heatmap_'+ca_ann[k]+'.png', bbox_inches='tight')
        plt.savefig(path_results_figures+option+'_heatmap_'+ca_ann[k]+'.pdf', bbox_inches='tight')
        if not opt_server:
            plt.show()
        plt.close()
        plt.clf()
        #print(np.min(CT_A_per_donor[k,:]))
        #print(np.max(CT_A_per_donor[k,:]))
        #print(np.mean(CT_A_per_donor[k,:]))

def plot_PCA_CT_percentages(DF_CT_A,ca_ann,path_results_figures):
    plt = loadPltSettings(20,20)
    for k in range(0,len(ca_ann)):
        data=DF_CT_A[DF_CT_A['method']==ca_ann[k]].pivot(index='Sample_ID',columns='CT_assigned', values='Annotation_percentage')
        pca = PCA(n_components=20)
        X=data.values
        pca.fit(X)
        print("Variance explained by all 20 principal components = ", sum(pca.explained_variance_ratio_*100))
        
        
        plt.plot(np.cumsum(pca.explained_variance_ratio_))
        plt.xlabel('Number of components')
        plt.ylabel('Explained variance')
        plt.savefig(path_results_figures+'PCA_explained_variance_'+ca_ann[k]+'.pdf', bbox_inches='tight')
        #plt.savefig(path_results_figures+'PCA_explained_variance_'+ca_ann[k]+'.png', bbox_inches='tight')
        plt.close()

        X_pca=pca.transform(X)
        group_vec = [np.unique(DF_CT_A[DF_CT_A['Sample_ID']==sample]['Group'])[0] for sample in data.index]
        PCA_R = pd.DataFrame(data = {'Sample_ID': data.index, 'PC1': X_pca[:,0], 'PC2': X_pca[:,1], 'Group': group_vec})

        plt.figure(figsize=(10,9))
        plt.rcParams['axes.labelsize']= 20
        fig, ax = plt.subplots()
        scatter_text(x_var='PC1',
                     y_var='PC2',
                     text_column='Sample_ID',
                     hue_variable='Group',
                     dataset=PCA_R,
                     title='PCA of cell type annotation percentages',
                     xlabel='PC 1 ('+str(np.round(pca.explained_variance_ratio_[0]*1000)/10)+' %)',
                     ylabel='PC 2 ('+str(np.round(pca.explained_variance_ratio_[1]*1000)/10)+' %)',
                     axes=ax,
                     opt_label_outliers_only = False)
        ax.set_legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,scatterpoints=1, fontsize=20)
        fig.savefig(path_results_figures+'PCA_'+ca_ann[k]+'.pdf', bbox_inches='tight')
        #plt.savefig(path_results_figures+'PCA_'+ca_ann[k]+'.png', bbox_inches='tight')
        plt.close(fig)

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

def plot_hist_max_contribution_per_cluster(df_max_per_cluster,path_results_figures,clustering_mode,opt_server):
    df_max_per_cluster['max contribution per cluster'] = df_max_per_cluster[0]*100
    plt = loadPltSettings(25,30)
    df_max_per_cluster['max contribution per cluster'].hist(bins=20,color='lightgrey')
    plt.xlabel('max contribution of a cell type per cluster')
    plt.ylabel('Absolute frequency')
    plt.ticklabel_format(axis="both", style="plain")
    plt.savefig(path_results_figures+'Max_cluster_contr_hist_'+clustering_mode+'.pdf', bbox_inches='tight')
    if not opt_server:
        plt.show()
    plt.close()
    plt.clf()

def plot_annotation_success_percentage(df,variables,path_results,opt_server):
    for var in variables:
        df1 = df.set_index(['n_genes','method']).copy()
        df0 = df1.reorder_levels(['n_genes','method']).sort_index()
        df0 = df0.unstack(level=-1) # unstack the 'Context' column
        fig, ax = plt.subplots()
        df0[var].plot(kind='bar', rot=0, ax=ax, color = {'scmap_cluster_HDO':'#8DADD7', 
                                                         'scmap_cluster_HVG_with_seurat_v3':  '#012F6C',
                                                         'scmap_cluster_HVG_with_cell_ranger': '#2762AF',
                                                         'scmap_cell2cluster_HDO':'#8DD7C5',
                                                         'scmap_cell2cluster_HVG_with_seurat_v3': '#006666',
                                                         'scmap_cell2cluster_HVG_with_cell_ranger':'#38907B'} )
        ax.legend(bbox_to_anchor=(1.1, 1.05), numpoints = 1)
        plt.ylabel(var)
        #plt.savefig(path_results+var+'_CT_annotations_per_method.png', bbox_inches='tight')
        plt.savefig(path_results+var+'_CT_annota_per_method.pdf', bbox_inches='tight')
        if not opt_server:
            plt.show()
        plt.close()
        plt.clf() 

def plot_UMAP_test_train(filename, path_data, path_results_figures):
    os.chdir(path_data)
    #plot UMAP of test and train set
    aD = anndata.read_loom(filename)
    # log-transformed and total count normalized data (each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization):
    sc.pp.normalize_total(aD, layers = 'all')
    sc.pp.log1p(aD,chunked=True) 
    os.chdir(path_results_figures)
    sc.pp.neighbors(aD, n_neighbors=10, n_pcs=30)
    sc.tl.umap(aD)
    v='group'
    sc.pl.umap(aD, color=v,save="_"+v+"_test_train")

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

def print_best_performing_CT_annotation_modes(DF_R):
    #print best performing methods
    unassigned_minimal = DF_R[DF_R['P_unassigned']==DF_R['P_unassigned'].min()][['method','n_genes']]
    correct_maximal = DF_R[DF_R['P_correct']==DF_R['P_correct'].max()][['method','n_genes']]
    false_minimal = DF_R[DF_R['P_false']==DF_R['P_false'].min()][['method','n_genes']]
    Cohens_kappa_maximal = DF_R[DF_R['Cohens_kappa']==DF_R['Cohens_kappa'].max()][['method','n_genes']]
    unassigned_CC_minimal = DF_R[DF_R['P_CC_unassigned']==DF_R['P_CC_unassigned'].min()][['method','n_genes']]
    correct_CC_maximal = DF_R[DF_R['P_CC_correct']==DF_R['P_CC_correct'].max()][['method','n_genes']]
    false_CC_minimal = DF_R[DF_R['P_CC_false']==DF_R['P_CC_false'].min()][['method','n_genes']]
    print("minimal % unassigned: \n" + str(unassigned_minimal)+'\n')
    print("maximal % correct: \n" + str(correct_maximal)+'\n')
    print("minimal % false: \n" + str(false_minimal)+'\n')
    print("maximal Cohen\'s kappa value: \n" + str(Cohens_kappa_maximal)+'\n')
    print("minimal % cell class unassigned: \n" + str(unassigned_CC_minimal)+'\n')
    print("maximal % cell class correct: \n" + str(correct_CC_maximal)+'\n')
    print("minimal % cell class false: \n" + str(false_CC_minimal)+'\n')

def print_difference_between_best_performing_methods(DF_R):
    #print difference between best performing methods
    d_kappa = DF_R[(DF_R['method']=='scmap_cell2cluster_HVG_with_cell_ranger') & (DF_R['n_genes']==5000)]['Cohens_kappa']-DF_R['Cohens_kappa'].max()
    d_unassigned = DF_R[(DF_R['method']=='scmap_cell2cluster_HVG_with_cell_ranger') & (DF_R['n_genes']==5000)]['P_unassigned']-DF_R['P_unassigned'].min()
    d_false = DF_R['P_false'].min()-DF_R[(DF_R['method']=='scmap_cluster_HVG_with_cell_ranger') & (DF_R['n_genes']==2000)]['P_false']
    d_correct = DF_R[(DF_R['method']=='scmap_cell2cluster_HVG_with_cell_ranger') & (DF_R['n_genes']==5000)]['P_correct']-DF_R['P_correct'].max()
    print("Difference between method scmap_cell2cluster_5000_HVG_with_cell_ranger and scmap_cluster_2000_HVG_with_cell_ranger:")
    print("difference in % unassigned: " + str(np.round(d_unassigned.iloc[0],2)))
    print("difference in % correct: " + str(np.round(d_correct.iloc[0],2)))
    print("difference in % false: " + str(np.round(d_false.iloc[0],2)))
    print("difference in Cohen\'s kappa value: " + str(np.round(d_kappa.iloc[0],2)))

def print_stats_volcano(DEGs_filtered,term):
    if "mito" in term:
        print("number of genes on pw "+term+":"+str(len(DEGs_filtered["Gene_short"].unique())))
        print("number of down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]<0))))
        print("Percentage of pw genes are down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]<0))/(len(DEGs_filtered["Gene_short"].unique())*len(DEGs_filtered["celltype"].unique()))))
        print("number of up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]>0))))
        print("Percentage of pw genes are up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]>0))/(len(DEGs_filtered["Gene_short"].unique())*len(DEGs_filtered["celltype"].unique()))))
        #are the MT genes part of the pathways?
        print("MT genes in term " + term + ":")
        print(DEGs_filtered[DEGs_filtered["Gene_short"].str.startswith("MT-")]["Gene_short"].unique())
    if term == "immune response":
        print("number of genes on pw "+term+":"+str(len(DEGs_filtered_micro)))
        print("number of down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered_micro["padj"]<=0.3,DEGs_filtered_micro["log2FoldChange"]<0))))
        print("Percentage of pw genes are down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered_micro["padj"]<=0.3,DEGs_filtered_micro["log2FoldChange"]<0))/len(DEGs_filtered_micro)))
        print("number of up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered_micro["padj"]<=0.3,DEGs_filtered_micro["log2FoldChange"]>0))))
        print("Percentage of pw genes are up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered_micro["padj"]<=0.3,DEGs_filtered_micro["log2FoldChange"]>0))/len(DEGs_filtered_micro)))
    if term=="neuron differentiation":
        print("number of genes on pw "+term+":"+str(len(DEGs_filtered["Gene_short"].unique())))
        print("number of down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]<0))))
        print("Percentage of pw genes are down DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]<0))/(len(DEGs_filtered["Gene_short"].unique())*len(DEGs_filtered["celltype"].unique()))))
        print("number of up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]>0))))
        print("Percentage of pw genes are up DEGs:"+str(np.sum(np.logical_and(DEGs_filtered["padj"]<=0.3,DEGs_filtered["log2FoldChange"]>0))/(len(DEGs_filtered["Gene_short"].unique())*len(DEGs_filtered["celltype"].unique()))))

