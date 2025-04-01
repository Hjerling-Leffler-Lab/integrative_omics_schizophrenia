
# -*- coding: utf-8 -*-
#Created on Wed Jun  3 14:38:41 2020
#@author: Lisa Bast
# version: 0.1

import numpy as np
import os
import loompy
import pandas as pd
import gc
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import seaborn as sns
from operator import itemgetter
import scrublet as scr
import anndata
import scanpy as sc
import math

FS=16
G=['CTRL','SCZ']
G2=['m','f']
add_str = '_cellranger'

## calculate:
def calc_stats_per_sample(D):
    sample_IDs,n_cells,disease,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
    min_CPG = np.empty(len(sample_IDs), dtype=float)
    mean_counts_per_cell = np.empty(len(sample_IDs), dtype=int)
    std_counts_per_cell = np.empty(len(sample_IDs), dtype=int)
    med_gpc = np.empty(len(sample_IDs), dtype=int)
    RLD = []
    for d_id,donor in enumerate(sample_IDs):
        mcg = calc_mean_counts_per_gene(D,donor)
        cpg = calc_cells_per_gene(D,donor)
        gene_list = (D.ra.Gene).astype('str').tolist()
        id_mito_genes = [i for i,item in enumerate(gene_list) if item.startswith('MT-')]
        id_ribo_genes = [i for i,item in enumerate(gene_list) if item.startswith('Rn4-')]
        genes_per_cell = calc_genes_per_cell(D,donor)
        counts_per_cell = calc_counts_per_cell(D,'',donor)
        spliced_counts_per_cell = []
        unspliced_counts_per_cell = []
        mito_counts_per_cell = calc_specific_counts_per_cell(id_mito_genes,D,donor)
        percent_mito = calc_percent_specific_counts(mito_counts_per_cell,counts_per_cell)
        ribo_counts_per_cell = calc_specific_counts_per_cell(id_ribo_genes,D,donor)
        percent_ribo = calc_percent_specific_counts(ribo_counts_per_cell,counts_per_cell)
        number_cells = len(D.ca.Donor[D.ca.Donor==donor])
        if d_id==0:
            MCG = [mcg]
            CPG = [cpg]
            GPC = [genes_per_cell]
            CPC = [counts_per_cell]
            sCPC = [spliced_counts_per_cell]
            usCPC = [unspliced_counts_per_cell]
            MCPC = [mito_counts_per_cell]
            PM = [percent_mito]
            RCPC = [ribo_counts_per_cell]
            PR = [percent_ribo]
            CPD = [number_cells]
        else:
            MCG.append(mcg)
            CPG.append(cpg)
            GPC.append(genes_per_cell)
            CPC.append(counts_per_cell)
            sCPC.append(spliced_counts_per_cell)
            usCPC.append(unspliced_counts_per_cell)
            MCPC.append(mito_counts_per_cell)
            PM.append(percent_mito)
            RCPC.append(ribo_counts_per_cell)
            PR.append(percent_ribo)
            CPD.append(number_cells)
        mean_counts_per_cell[d_id] = np.mean(CPC[d_id])
        std_counts_per_cell[d_id] = np.std(CPC[d_id])
        min_CPG[d_id] = np.min(CPG[d_id])  
        med_gpc[d_id] = np.median(GPC[d_id])
    return np.array(MCG), np.array(CPG), np.array(min_CPG), np.array(GPC), np.array(CPC), np.array(sCPC), np.array(usCPC), np.array(MCPC), np.array(PM), np.array(RCPC), np.array(PR), np.array(mean_counts_per_cell), np.array(std_counts_per_cell), np.array(RLD), np.array(med_gpc), np.array(CPD)

def calc_mean_counts_per_gene(D,donor):
    idx = np.where(D.ca.Donor==donor)
    mcg = np.mean(D[:,idx[0]],axis=1)
    return mcg

def calc_percent_specific_counts(specific_counts_per_cell,counts_per_cell):
    percent_specific_counts = 100*np.divide(specific_counts_per_cell,counts_per_cell)
    return percent_specific_counts

def calc_cells_per_gene(D,donor):
    if donor=='any':
        sample_IDs = np.unique(D.ca.Donor)
        for s_id,donor in enumerate(sample_IDs):
            idx = np.where(D.ca.Donor==donor)
            if s_id==0:
                cpg = np.sum(D[:,idx[0]]>0,axis=1)
            else:
                cpg = cpg + np.sum(D[:,idx[0]]>0,axis=1)
        #cpg = np.sum(D[:,:],axis=1)
    else:
        idx = np.where(D.ca.Donor==donor)
        cpg = np.sum(D[:,idx[0]]>0,axis=1)
    return cpg

def calc_genes_per_cell(D,donor):
    idx = np.where(D.ca.Donor==donor)
    genes_per_cell = np.sum(D[:,idx[0]]>0,axis=0)
    return genes_per_cell

def calc_counts_per_cell(D,layer,donor):
    idx = np.where(D.ca.Donor==donor)
    counts_per_cell = np.sum(D[layer][:,idx[0]],axis=0)
    return counts_per_cell

def calc_specific_counts_per_cell(id_specific_genes,D,donor):
    idx = np.where(D.ca.Donor==donor)
    specific_counts_per_cell = np.sum(itemgetter(id_specific_genes)(D[:,idx[0]]),axis=0)
    return specific_counts_per_cell

def calc_doublet_prediction_scrublet(donor, add_str, add_str_pagoda, TH_min_number_of_cells_expr_a_gene, TH_min_number_counts, path_data_filt, path_results, opt_plot):
    gc.collect()
    os.chdir(path_data_filt)
    with loompy.connect(donor+add_str+add_str_pagoda+'_cell_TH_filtered.loom') as D:
        counts_matrix = np.asmatrix(D[:,:].transpose())
        print(np.shape(counts_matrix))
    #np.shape(counts_matrix)
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate = 0.05)
    #TO DO: make sure to filter out samples with too few cells
    #otherwise ValueError: n_components=30 must be between 0 and min(n_samples, n_features)=18 with svd_solver='full'
    counts_matrix = None
    doublet_scores, predicted_doublets = scrub.scrub_doublets()#min_counts=TH_min_number_counts, min_cells=TH_min_number_of_cells_expr_a_gene, min_gene_variability_pctl=85, n_prin_comps=30)
    if opt_plot==True:
        plot_doublets_scrublet(scrub,donor,path_results)
    #doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=TH_min_counts[d_id], min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    #print(doublet_scores)
    #print(predicted_doublets)
    return scrub, doublet_scores, predicted_doublets 

def calculate_spearman_corr_with_QC_metrics(df_qc_and_pca):
    #spearman correlation for PCx with QCy
    qc_ids = range([i for i in range(len(df_qc_and_pca.columns)) if df_qc_and_pca.columns[i] == 'num_umi'][0],len(df_qc_and_pca.columns))
    #pc_ids = np.where(df_qc_and_pca.columns.str.startswith('PC'))
    pc_ids = range(np.min(np.where(df_qc_and_pca.columns.str.startswith('PC'))),np.max(np.where(df_qc_and_pca.columns.str.startswith('PC'))))
    Corr_spearman = np.zeros((np.max(np.shape(qc_ids)),np.max(np.shape(pc_ids))))
    #pval_spearman = np.empty((np.max(np.shape(qc_ids)),np.max(np.shape(pc_ids))), dtype=float)
    for j,pc_id in enumerate(pc_ids):
        for i,qc_id in enumerate(qc_ids):
            a=np.array(df_qc_and_pca[df_qc_and_pca.columns[qc_ids][i]].rank(), dtype ='float')
            b=np.array(df_qc_and_pca[df_qc_and_pca.columns[pc_ids][j]].rank(), dtype ='float')
            corr_ab = np.corrcoef(x=a,y=b)[0,1]
            if corr_ab<1 and corr_ab>-1:
                Corr_spearman[i,j] = np.around(corr_ab,4)
            else:
                Corr_spearman[i,j] = math.nan    
    return Corr_spearman, qc_ids

## create:
def create_loom_file(D, genes_kept, cells_kept, path,filename):
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

def create_TH_filtered_loom_files(sample_IDs, disease_status,  opt_pagoda_filtering, path_data_raw,path_data_filt,path_code,path_results):
    if opt_pagoda_filtering:
        add_str_pagoda = "_pagoda"
        add_str_pagoda_long = "_pagoda_filtered"
        path_data_source = path_data_filt
    else:
        add_str_pagoda = ""
        add_str_pagoda_long = ""
        path_data_source = path_data_raw
    N_cells_discarded_1 = np.empty((len(sample_IDs),1), dtype=int)
    N_cells_discarded_2 = np.empty((len(sample_IDs),1), dtype=int)
    N_cells_discarded_3 = np.empty((len(sample_IDs),1), dtype=int)
    N_cells_discarded_4 = np.empty((len(sample_IDs),1), dtype=int)
    N_cells_discarded_5 = np.empty((len(sample_IDs),1), dtype=int)
    N_cells_discarded_total = np.empty((len(sample_IDs),1), dtype=int)
    #N_genes_discarded_per_sample = np.empty((len(sample_IDs),1), dtype=int)
    #To do:
    #N_cells_discarded_5 = np.empty((len(sample_IDs),1), dtype=int)
    P_genes_discarded = np.empty((len(sample_IDs),1), dtype=float)
    P_cells_discarded = np.empty((len(sample_IDs),1), dtype=float)
    sample_IDs_after_TH_filtering = list()
    sample_IDs_removed_by_TH_filtering = list()
    bool_first = True
    TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_genes, TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_max_ribosomal_fraction, TH_min_cells_per_sample = get_TH_values()
    
    os.chdir(path_data_filt)
    source_filename = "Samples"+add_str+add_str_pagoda_long+".loom"
    
    with loompy.connect(source_filename) as D_all:
        N_genes_discarded_overall, genes_to_remove_overall = get_rare_genes(D_all,'any',TH_min_number_of_cells_expr_a_gene_all_samples)
        genes_to_remove_overall_str = D_all.ra['Gene'][tuple(genes_to_remove_overall)]
        genes_kept = genes_to_remove_overall[0]==False
    for d_id,donor in enumerate(sample_IDs):
        os.chdir(path_data_source)
        with loompy.connect(donor+add_str+add_str_pagoda_long+".loom") as D:
            genes_kept = [all(sg != gr for gr in genes_to_remove_overall_str) for sg in D.ra['Gene']]
            cells_kept = np.ones((1, np.shape(D)[1]), dtype=bool).flatten()
            P_genes_discarded[d_id] = 100*np.sum(genes_kept==False)/len(genes_kept)
            filename = donor+add_str+add_str_pagoda+'_gene_TH_filtered.loom'
            if np.any(genes_kept) and np.sum(genes_kept)>30:
                create_loom_file(D, genes_kept, cells_kept,path_data_filt,filename)
                print('Sample ' + donor + ' is gene threshold filtered.')
            else:
                if not np.any(genes_kept):
                    reason = ' all genes have been removed!'
                elif np.sum(genes_kept)<30:
                    reason = ' too few genes are left!'
                print('Sample ' + donor + ' is removed during threshold filtering as' + reason)
                sample_IDs_removed_by_TH_filtering.append(donor)
    #find cells to discard based on single samples:
    target_files_CTRL=[]
    target_files_SCZ=[]
    for d_id,donor in enumerate(sample_IDs):
        os.chdir(path_data_filt) 
        with loompy.connect(donor+add_str+add_str_pagoda+'_gene_TH_filtered.loom') as D:
            #cell thresholds per sample:
            N_cells_discarded_1[d_id], cells_kept_1 = get_cells_with_many_genes(D,donor,TH_min_number_genes)
            N_cells_discarded_2[d_id], cells_kept_2 = get_cells_with_low_mito_fraction(D,donor,TH_max_fraction_mito_counts)
            N_cells_discarded_3[d_id], cells_kept_3 = get_cells_with_low_ribosomal_fraction(D,donor,TH_max_ribosomal_fraction)
            N_cells_discarded_4[d_id], cells_kept_4 = get_cells_with_high_count_depth(D,donor,"",TH_min_number_counts)
            N_cells_discarded_total[d_id] = N_cells_discarded_1[d_id]+N_cells_discarded_2[d_id]+N_cells_discarded_3[d_id]+N_cells_discarded_4[d_id]
            cells_kept = np.logical_and(np.logical_and(np.logical_and(cells_kept_1[0],cells_kept_2[0]),cells_kept_3[0]),cells_kept_4[0])
            
            P_cells_discarded[d_id] = 100*(np.sum(cells_kept==False)/len(cells_kept))#100*(N_cells_discarded_1[d_id]+N_cells_discarded_2[d_id]+N_cells_discarded_3[d_id]+N_cells_discarded_4[d_id])/np.shape(D)[1]
            if np.any(cells_kept) and np.sum(cells_kept)>TH_min_cells_per_sample:
                sample_IDs_after_TH_filtering.append(donor) 
                filename = donor+add_str+add_str_pagoda+'_cell_TH_filtered.loom'
                genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
                create_loom_file(D, genes_kept, cells_kept, path_data_filt,filename)
                if bool_first:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = [filename]
                    else:
                        target_files_SCZ = [filename]
                    bool_first = False
                else:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = target_files_CTRL + [filename]
                    else:
                        target_files_SCZ = target_files_SCZ + [filename]
                print('Sample ' + donor + ' is cell threshold filtered!')
            else:
                if not np.any(cells_kept):
                    reason = ' all cells have been removed!'
                elif np.sum(cells_kept)<TH_min_cells_per_sample:
                    reason = ' too few cells are left!'
                print('Sample ' + donor + ' is removed during threshold filtering as' + reason)
                sample_IDs_removed_by_TH_filtering.append(donor)
    ### Data frame with statistics:
    #, 'Mean Counts per cell after filtering':
    df_filtering_result = pd.DataFrame(data={'Donor ID': sample_IDs, 
                                             'Disease': disease_status,
                                             'Number of cells removed total': N_cells_discarded_total.flatten(),
                                             'Cells <='+str(TH_min_number_genes)+' genes': N_cells_discarded_1.flatten(), 
                                             'Cells <='+str(TH_min_number_counts)+' counts': N_cells_discarded_4.flatten(), 
                                             'Cells >'+str(TH_max_fraction_mito_counts*100)+'% mito counts': N_cells_discarded_2.flatten(), 
                                             'Cells >'+str(TH_max_ribosomal_fraction*100)+'% ribo counts': N_cells_discarded_3.flatten(),
                                             #'Number of genes removed total': N_genes_discarded_per_sample.flatten(), 
                                             'Genes overall expressed in <'+str(TH_min_number_of_cells_expr_a_gene_all_samples)+ ' cells':N_genes_discarded_overall*np.ones((len(sample_IDs),1)).flatten(),
                                             #'Genes expressed in <'+str(TH_min_number_of_cells_expr_a_gene)+ ' cells':N_genes_discarded_1.flatten(),
                                             'Percentage Genes discarded': P_genes_discarded.flatten(), 
                                             'Percentage Cells discarded': P_cells_discarded.flatten()})
    #create path if does not already exist:
    isExist = os.path.exists(path_results)
    if not isExist:    
        os.makedirs(path_results)
    os.chdir(path_results)
    df_filtering_result.to_pickle('df'+add_str_pagoda+'_TH_filtering_result'+add_str)
    df_filtering_result.to_excel('df'+add_str_pagoda+'_TH_filtering_result'+add_str+'.xlsx')
    os.chdir(path_data_filt)
    ### create loom file with all/ all CTRL/ all SCZ samples
    if target_files_CTRL:
        print('target_files CTRL:')
        print(target_files_CTRL)
        loompy.combine(target_files_CTRL, 'Samples_CTRL'+add_str+add_str_pagoda+'_TH_filtered.loom', batch_size = 3000)
    if target_files_SCZ:
        print('target_files SCZ:')
        print(target_files_SCZ)
        loompy.combine(target_files_SCZ, 'Samples_SCZ'+add_str+add_str_pagoda+'_TH_filtered.loom', batch_size = 3000)
    loompy.combine(files = target_files_CTRL+target_files_SCZ, output_file = 'Samples'+add_str+add_str_pagoda+'_TH_filtered.loom', batch_size = 3000)   
    os.chdir(path_results)
    ### save which samples are removed
    np.save('sample_IDs_after_TH_filtering'+add_str+add_str_pagoda,sample_IDs_after_TH_filtering)
    np.save('sample_IDs_removed_by_TH_filtering'+add_str+add_str_pagoda,sample_IDs_removed_by_TH_filtering)
    gc.collect()
    os.chdir(path_code)
    return df_filtering_result

def create_and_merge_loom_files(path_data,path_patient_info_data,path_summary_metrics,path_code,library_nrs,samples_to_drop):

    folder = ""
        
    os.chdir(path_data)
    SI = pd.read_excel(path_patient_info_data+"T1_basic_donor_information.xlsx",'basic_donor_information')
    #filter for snRNAseq samples:
    SI = SI[SI["snRNAseq_performed"].str.startswith("yes")]
    counter = 1
    for l_id,lib in enumerate(library_nrs):
        lib_sel_bools = SI["donor_ID_internal_7"].str.startswith('10x'+str(library_nrs[l_id])+'_')
        lib_sample_nr = SI["donor_ID_internal_7"][lib_sel_bools==True]
        sample_ID_python = SI["donor_ID_internal_8"][lib_sel_bools==True]
        groups = SI["scz_status"][lib_sel_bools==True]
        groups[groups=="scz"]="SCZ"
        groups[groups=="ctrl"]="CTRL"
        for s2d in samples_to_drop:
            #drop sample 10x2_19 and others as specified
            s2d_list = lib_sample_nr[lib_sample_nr==s2d].index.values.tolist()
            if len(s2d_list)!=0:
                lib_sample_nr = lib_sample_nr.drop(s2d_list)
                sample_ID_python = sample_ID_python.drop(s2d_list)
                groups  = groups.drop(s2d_list)
        sample_folder = ["./library_"+str(lib)+folder+"/Counts_"+str(lib_sample_nr.iloc[i])+".loom" for i in range(0,len(lib_sample_nr))] 
        
        #print(lib_sample_nr)
        for i in range(0,len(lib_sample_nr)):
            if l_id==0 and i==0:
                target_files = [sample_ID_python.iloc[i]+add_str+".loom"]
            else:
                target_files = target_files + [sample_ID_python.iloc[i]+add_str+".loom"]
            load_aligned_data_and_metadata(sample_folder[i],sample_ID_python.iloc[i]+add_str,groups.iloc[i],'10x'+str(lib),path_data,path_patient_info_data,path_summary_metrics)
            counter=counter+1
    os.chdir(path_data)
    print(target_files)
    loompy.combine(files = target_files, output_file = 'Samples'+add_str+'.loom', key= 'Gene', batch_size = 3000)
    os.chdir(path_code)

def create_TH_and_D_filtered_loom_files(sample_IDs,disease_status,PD,opt_manual_doublet_TH_adjustment,opt_pagoda_filtering,add_str,TH_min_cells_per_sample,path_data_raw, path_data_filt,path_code,path_results,opt_D_version):
    os.chdir(path_data_filt)
    bool_first = True
    sample_IDs_after_D_filtering = list()
    sample_IDs_removed_by_D_filtering = list()
    #disease = list()
    N_cells_discarded = np.empty((len(sample_IDs),1), dtype=int)
    P_cells_discarded = np.empty((len(sample_IDs),1), dtype=float)
    target_files_CTRL=[]
    target_files_SCZ=[]
    if opt_manual_doublet_TH_adjustment:
        add_str2 = '_adj'+'_'+opt_D_version
    else:
        add_str2 = ''
        
    if opt_pagoda_filtering:
        add_str_pagoda = "_pagoda"
    else:
        add_str_pagoda = ""
        
    for d_id,donor in enumerate(sample_IDs):
        with loompy.connect(donor+add_str+add_str_pagoda+'_cell_TH_filtered.loom') as D_TH_filtered:
            #disease.append(np.unique(D_TH_filtered.ca.Disease[np.where(D_TH_filtered.ca.Donor==donor)])[0])
            #PD is a list containing np arrays of bools for every donor indicating if cell is predicted to be a doublet or not
            N_cells_discarded[d_id], cells_kept = get_predicted_non_doublets(PD[d_id])
            P_cells_discarded[d_id] = (100*N_cells_discarded[d_id])/np.shape(D_TH_filtered)[1]
            print('number of cells in sample' + donor + ': ' + str(np.sum(cells_kept)))
            if np.sum(cells_kept)>=TH_min_cells_per_sample:#np.any(cells_kept):
                filename=donor+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom'
                genes_kept = np.ones((1, np.shape(D_TH_filtered)[0]), dtype=bool).flatten()
                create_loom_file(D_TH_filtered, genes_kept, cells_kept, path_data_filt,filename)
                sample_IDs_after_D_filtering.append(donor)
                if bool_first:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = [donor+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom']
                    else:
                        target_files_SCZ = [donor+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom']
                    bool_first = False
                else:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = target_files_CTRL + [donor+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom']
                    else:
                        target_files_SCZ = target_files_SCZ + [donor+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom']
                print('Sample ' + donor + ' is threshold and doublet filtered!')
            else:
                print('Sample ' + donor + ' is removed during doublet filtering!')
                sample_IDs_removed_by_D_filtering.append(donor)
    ### create loom file with all/ all CTRL/ all SCZ samples
    if target_files_CTRL:
        print('target_files_CTRL:')
        print(target_files_CTRL)
        loompy.combine(files = target_files_CTRL, output_file = 'Samples_CTRL'+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom', batch_size = 3000)
    if target_files_SCZ:
        print('target_files_SCZ:')
        print(target_files_SCZ)
        loompy.combine(files = target_files_SCZ, output_file = 'Samples_SCZ'+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom', batch_size = 3000)
    loompy.combine(files = target_files_CTRL+target_files_SCZ, output_file = 'Samples'+add_str+add_str_pagoda+'_TH_and_D'+add_str2+'_filtered.loom', batch_size = 3000)   
    os.chdir(path_results)
    ### save which samples are removed
    np.save('sample_IDs_after_D_filtering'+add_str+add_str_pagoda,sample_IDs_after_D_filtering)
    np.save('sample_IDs_removed_by_D_filtering'+add_str+add_str_pagoda,sample_IDs_removed_by_D_filtering)
    ### Data frame with statistics:
    df_filtering_result = pd.DataFrame(data={'Donor ID': sample_IDs, 'Disease': disease_status, 'Number of cells removed (predicted doublets)': N_cells_discarded.flatten(),  'Percentage Cells discarded': P_cells_discarded.flatten()})
    df_filtering_result.to_pickle('df'+add_str_pagoda+'_TH_and_D'+add_str2+'filtering_result'+add_str)
    df_filtering_result.to_excel('df'+add_str_pagoda+'_TH_and_D'+add_str2+'filtering_result'+add_str+'.xlsx')
    os.chdir(path_code)
    return df_filtering_result

def integrate_filtering_status_in_metadata(path_data_raw,path_data_filt, opt_D_version, opt_pagoda_filtering):
    if opt_pagoda_filtering:
        add_str_pagoda = '_pagoda'
    else:
        add_str_pagoda = ''
    #get sample IDs of samples post filtering:    
    filename_raw = 'Samples'+add_str+'.loom'
    filename_b = 'Samples' + add_str + add_str_pagoda + '_TH_filtered.loom'
    filename_c = 'Samples' + add_str + add_str_pagoda + '_TH_and_D_adj_' + opt_D_version + '_filtered.loom'#_and_CT_annotated.loom'
    if opt_pagoda_filtering == True:
        filename_a = 'Samples'+add_str+'_pagoda_filtered.loom'
        with loompy.connect(path_data_filt+filename_a) as D_a:
            samples_a,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D_a)
    with loompy.connect(path_data_filt+filename_b) as D_b:
        samples_b,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D_b) 
    with loompy.connect(path_data_raw+filename_raw) as D_raw:  
        samples_raw,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D_raw)  
    with loompy.connect(path_data_filt+filename_c) as D_c:
        samples_c,_,_,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D_c)
    
    if opt_pagoda_filtering == True:
        samples_after_filtering = np.intersect1d(np.intersect1d(samples_a,samples_b),samples_c)
    else:
        samples_after_filtering = np.intersect1d(samples_b,samples_c)
        
    #for every sample in the raw data:
    #check for each sample which cells got filtered --> add to status variable in col metadata
    for d_id,donor in enumerate(samples_raw):
        filename_raw_i = donor+add_str+'.loom'
        filename_b_i = donor+add_str+add_str_pagoda+'_cell_TH_filtered.loom'
        filename_c_i = donor+add_str+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered.loom'#'_and_CT_annotated.loom'
        D_raw_i = loompy.connect(path_data_raw+filename_raw_i)
        all_cellIDs = D_raw_i.ca.CellID.flatten()
        if donor in samples_after_filtering:
            new_col_attr = np.repeat(['not filtered'],len(all_cellIDs))#np.shape(D_b[:,:])[1])
            #for TH and doublet filtering get difference:       
            D_b_i = loompy.connect(path_data_filt+filename_b_i)#TH filtered
            cellIDs_TH_filtered = D_b_i.ca.CellID.flatten()
            D_b_i.close()
            D_c_i = loompy.connect(path_data_filt+filename_c_i)#Doublet filtered
            cellIDs_TH_and_D_filtered = D_c_i.ca.CellID.flatten()
            D_c_i.close()
            TH_filered_cells = np.setdiff1d(all_cellIDs,cellIDs_TH_filtered,assume_unique=True)
            D_filtered_cells = np.setdiff1d(cellIDs_TH_filtered,cellIDs_TH_and_D_filtered,assume_unique=True)
            #go only through those items
            for cid,c_id in enumerate(TH_filered_cells):
                new_col_attr[all_cellIDs==c_id]='TH filtered'
            for cid,c_id in enumerate(D_filtered_cells):
                new_col_attr[all_cellIDs==c_id] = 'D filtered'
                        #in case pagoda filtering was performed get difference:
            if opt_pagoda_filtering == True:
                filename_a_i = donor+add_str+add_str_pagoda+'_filtered.loom'
                D_a_i = loompy.connect(path_data_filt+filename_a_i)
                cellIDs_pagoda_filtered = D_a_i.ca.CellID.flatten()
                D_a_i.close()
                Pagoda_filtered_cells = np.setdiff1d(all_cellIDs,cellIDs_pagoda_filtered,assume_unique=True)
                for cid,c_id in enumerate(Pagoda_filtered_cells):
                    new_col_attr[all_cellIDs==c_id] = 'P filtered'
        else: 
            new_col_attr = np.repeat(['Sample removed'],len(all_cellIDs))
        new_filename_i = add_col_metadata_to_loom_file(D_raw_i, new_col_attr, 'filtering_status', path_data_raw, filename_raw_i[:-5]+add_str_pagoda+'.loom', True)
        if d_id==0:
            new_filenames = [new_filename_i]
        else:
            new_filenames.append(new_filename_i)
    loompy.combine(files = new_filenames, output_file = 'Samples'+add_str+'.loom', batch_size = 3000)
    return filename_raw

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

## get:
def get_pagoda_filtered_data(sample_IDs, disease_status, path_code, path_data_filt, path_results, opt_create_and_merge_raw_data_loom_files, opt_use_TH_filtered_data, opt_use_TH_and_D_filtered_data, opt_plot_filtered_data, TH_min_cells_per_sample,opt_D_version, opt_use_pagoda_filtered_data):

    if opt_use_pagoda_filtered_data:
        os.chdir(path_results)
        ### Data frame with statistics:
        df_filtering_result= pd.read_pickle('df_pagoda_filtering_result'+add_str)
        os.chdir(path_data_filt)
    else:
        bool_first = True
        sample_IDs_after_pagoda_filtering = list()
        sample_IDs_removed_by_pagoda_filtering = list()
        #disease = list()
        N_cells_discarded = np.empty((len(sample_IDs),1), dtype=int)
        P_cells_discarded = np.empty((len(sample_IDs),1), dtype=float)
        target_files_CTRL=[]
        target_files_SCZ=[]
        for d_id,donor in enumerate(sample_IDs):
            #load raw data
            os.chdir(path_data_filt)
            S = loompy.connect(donor+add_str+".loom")
            #load filtered data
            os.chdir(path_data_filt)
            S_pf = loompy.connect(donor+add_str+"_pagoda_filtered_R.loom",validate=False)
            #store number of genes and cells
            n_genes_S, n_cells_S = np.shape(S)
            #which cells are kept?
            cells_kept = S_pf.ca["Pagoda filtered"]==0
            #check if CellIDs match
            if not(np.sum(S_pf.ca["CellID"]==S.ca["CellID"])==n_cells_S):
                print("files "+donor+add_str+".loom"+" and "+ donor+add_str+"_pagoda_filtered_R.loom"+" contain different cells/ a different order.")
            #check if genes match
            if not(np.sum(S.ra["Gene"]==S_pf.ra["Gene"])==n_genes_S):
                print("files "+donor+add_str+".loom"+" and "+ donor+add_str+"_pagoda_filtered_R.loom"+" contain different genes/ a different order.")
            #identify cells present in both files
            N_cells_discarded[d_id]= np.sum(S_pf.ca['Pagoda filtered'])
            S_pf.close()
            P_cells_discarded[d_id] = (100*N_cells_discarded[d_id])/n_cells_S
            print('number of cells in sample ' + donor + ': ' + str(np.sum(cells_kept)))
            if np.sum(cells_kept)>=TH_min_cells_per_sample:#np.any(cells_kept):
                filename=donor+add_str+'_pagoda_filtered.loom'
                genes_kept = np.ones((1, np.shape(S)[0]), dtype=bool).flatten()
                create_loom_file(S, genes_kept, cells_kept, path_data_filt,filename)
                S.close()
                sample_IDs_after_pagoda_filtering.append(donor)
                if bool_first:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = [donor+add_str+'_pagoda_filtered.loom']
                    else:
                        target_files_SCZ = [donor+add_str+'_pagoda_filtered.loom']
                    bool_first = False
                else:
                    if disease_status[d_id]=='CTRL':
                        target_files_CTRL = target_files_CTRL + [donor+add_str+'_pagoda_filtered.loom']
                    else:
                        target_files_SCZ = target_files_SCZ + [donor+add_str+'_pagoda_filtered.loom']
                print('Sample ' + donor + ' is pagoda filtered!')
            else:
                print('Sample ' + donor + ' is removed during pagoda filtering!')
                sample_IDs_removed_by_pagoda_filtering.append(donor)
        ### create loom file with all/ all CTRL/ all SCZ samples
        if target_files_CTRL:
            print('target_files_CTRL:')
            print(target_files_CTRL)
            loompy.combine(files = target_files_CTRL, output_file = 'Samples_CTRL'+add_str+'_pagoda_filtered.loom', batch_size = 3000)
        if target_files_SCZ:
            print('target_files_SCZ:')
            print(target_files_SCZ)
            loompy.combine(files = target_files_SCZ, output_file = 'Samples_SCZ'+add_str+'_pagoda_filtered.loom', batch_size = 3000)
        loompy.combine(files = target_files_CTRL+target_files_SCZ, output_file = 'Samples'+add_str+'_pagoda_filtered.loom', batch_size = 3000)  
             
        os.chdir(path_results)
        ### save which samples are removed
        np.save('sample_IDs_after_pagoda_filtering'+add_str,sample_IDs_after_pagoda_filtering)
        np.save('sample_IDs_removed_by_pagoda_filtering'+add_str,sample_IDs_removed_by_pagoda_filtering)
        ### Data frame with statistics:
        df_filtering_result = pd.DataFrame(data={'Donor ID': sample_IDs, 'Disease': disease_status, 'Number of cells removed (pagoda filtering)': N_cells_discarded.flatten(),  'Percentage Cells discarded': P_cells_discarded.flatten()})
        df_filtering_result.to_pickle('df_pagoda_filtering_result'+add_str)
        df_filtering_result.to_excel('df_pagoda_filtering_result'+add_str+'.xlsx')

    os.chdir(path_data_filt)
    #load pagoda_filtered loom files
    D_fil_CTRL = loompy.connect("Samples_CTRL"+add_str+"_pagoda_filtered.loom")
    D_fil_SCZ = loompy.connect("Samples_SCZ"+add_str+"_pagoda_filtered.loom")
    D_fil = loompy.connect("Samples"+add_str+"_pagoda_filtered.loom")
    os.chdir(path_code)
    return D_fil, D_fil_CTRL, D_fil_SCZ, df_filtering_result

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

def get_TH_values():
    #TH_min_cell_cluster_size_percent = 0.025 #smalles cell cluster has 2,5% of cells in the sample assigned
    TH_min_number_of_cells_expr_a_gene_all_samples = 1 #based on cells of all samples
    #TH_min_number_of_cells_expr_a_gene = 1 #based on cells of 1 sample
    #TH_min_number_of_s_counts_per_gene = [] #not taken into account right now
    #TH_min_number_of_us_counts_per_gene = [] #not taken into account right now
    TH_min_number_genes = 500 #expressed in a cell
    TH_min_number_counts = 2000 # TH for microglia, for neurons: set to 10000 later; with pagoda 500
    TH_min_number_s_counts = []
    TH_min_number_us_counts = []
    TH_max_fraction_mito_counts= 1
    #To do:
    #TH_min_protein_coding_fraction = 0
    TH_max_ribosomal_fraction = 1
    TH_min_cells_per_sample = 500 # after applying thresholds this is the number of cells that a sample must contain in order to get further analysed
    #return TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_of_cells_expr_a_gene, TH_min_number_genes, TH_min_number_counts, TH_min_number_of_s_counts_per_gene, TH_min_number_of_us_counts_per_gene, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_max_ribosomal_fraction, TH_min_cells_per_sample
    return TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_genes, TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_max_ribosomal_fraction, TH_min_cells_per_sample

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


def get_paths(path_project,opt_pagoda_filtering,opt_results_path,n_cluster):
    if opt_results_path=='raw':
        path_data = path_project + "/3_quality_control/output/raw/"
    elif opt_results_path=='preprocessed':
        path_data = path_project + "/3_quality_control/output/filtered/"
    path_results = path_data
    path_code = path_project + "/3_quality_control/script/"
        
    return path_data, path_code, path_results

def get_predicted_non_doublets(PD_bools):
    cells_kept = PD_bools==False
    n_cells_discarded = len(PD_bools)-np.sum(cells_kept)
    return n_cells_discarded, cells_kept 

def get_color_maps_for_anndata(aD):
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_76CTs_colors"] = 'nipy_spectral'
    aD.uns["CT_ann_ABM_MCA_scmap_cluster_51_CTs_colors"] = 'nipy_spectral'
    aD.uns["Disease_colors"] = ['#487BAF','#FA9835']#['steelblue','darkorange']
    aD.uns["Sex_colors"] = ['#CB5C93','#0080FF']#['palevioletred','dodgerblue']
    aD.uns["filtering_status_colors"] = ["#CC0066","#078383","#939393","#2284C6","#000000"]
    return aD

def load_aligned_data_and_metadata(pathfile,fileID,disease,library_number,path_data,path_patient_info_data,path_summary_metrics):
    # Created on Tue Apr 21 15:59:36 2020
    # @author: José Martinez Lopez, adapted by Lisa Bast
    with loompy.connect(pathfile) as s:
        total = s[:,:]
    
        dt = np.dtype('U20')
        cellid=np.array(np.repeat(fileID[0:fileID.find('_')],repeats=total.shape[1]),dtype=dt)
        
        for i in range(0,total.shape[1]):
            cellid[i] =str(i)+cellid[i]
        
        filename =fileID+'.loom'
        
        genenames = s.ra["Gene"]
        accession = s.ra["Accession"]

        row_attrs = {"Gene": genenames,
                    "Accession": accession}
    donor = fileID[0:fileID.find('_')]
    print(donor)
    
    current_path = os.getcwd()
    donor_info = pd.read_excel(path_patient_info_data+'sample_info_SUN_SCZ.xlsx', 'sample_info') 
    qc_metrics = pd.read_csv(path_summary_metrics+'metrics_summary_tidy.csv')
    os.chdir(current_path)
    sex = donor_info[donor_info['donor_ID_internal_8']==donor]['sex'].iloc[0]
    age = donor_info[donor_info['donor_ID_internal_8']==donor]['age'].iloc[0]
    pmi_h = donor_info[donor_info['donor_ID_internal_8']==donor]['PMI_h'].iloc[0]
    if type(pmi_h) == int:
        pmi_h = float(pmi_h)
    elif type(pmi_h) == str:
        pmi_h = float(pmi_h[0:2])+float(pmi_h[3:5])/60+float(pmi_h[6:8])/(60*60)
    elif type(pmi_h) != float and type(pmi_h) != str and type(pmi_h) != int:
        print('type of pmi_h is:' + str(type(pmi_h)))
    if type(age) == str and age.strip()=='':
        age = float("nan")
    else:
        age=float(age)
    if len(qc_metrics[qc_metrics['donor_ID_python']==donor]['num_umi'])>0:
        num_umi = qc_metrics[qc_metrics['donor_ID_python']==donor]['num_umi'].iloc[0]
        print("number umis: "+str(num_umi))
        mean_reads_per_umi = qc_metrics[qc_metrics['donor_ID_python']==donor]['mean_reads_per_umi'].iloc[0]
        num_reads = qc_metrics[qc_metrics['donor_ID_python']==donor]['num_reads'].iloc[0]
        p_valid_barcodes = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_valid_barcodes'].iloc[0]
        p_sequencing_saturation = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_sequencing_saturation'].iloc[0]
        p_genome_not_gene = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_genome_not_gene'].iloc[0]
        p_mapped_reads = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_mapped_reads'].iloc[0]
        p_unmapped_reads = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_unmapped_reads'].iloc[0]
        if 'p_cell_del_filt' in qc_metrics.columns:
            p_cell_del_filt = qc_metrics[qc_metrics['donor_ID_python']==donor]['p_cell_del_filt'].iloc[0]
        else:
            p_cell_del_filt = float('NaN')
        if 'mean_counts_per_barcode' in qc_metrics.columns:
            mean_counts_per_barcode = qc_metrics[qc_metrics['donor_ID_python']==donor]['mean_counts_per_barcode'].iloc[0]
        else:
            mean_counts_per_barcode = float('NaN')
        if 'std_counts_per_barcode' in qc_metrics.columns:
            std_counts_per_barcode = qc_metrics[qc_metrics['donor_ID_python']==donor]['std_counts_per_barcode'].iloc[0]
        else:
            std_counts_per_barcode = float('NaN')
        if 'median_gpc' in qc_metrics.columns:
            median_gpc = qc_metrics[qc_metrics['donor_ID_python']==donor]['median_gpc'].iloc[0]
        else:
            median_gpc = float('NaN')
    else:
        num_umi = float('NaN')
        mean_reads_per_umi = float('NaN')
        num_reads = float('NaN')
        p_valid_barcodes = float('NaN')
        p_sequencing_saturation = float('NaN')
        p_genome_not_gene = float('NaN')
        p_mapped_reads = float('NaN')
        p_unmapped_reads = float('NaN')
        p_cell_del_filt = float('NaN')
        mean_counts_per_barcode = float('NaN')
        std_counts_per_barcode = float('NaN')
        median_gpc = float('NaN')
    col_attrs = { "CellID": cellid, 
                  "Disease": np.repeat(disease,repeats=total.shape[1]), #assign disease status
                  "Donor": np.repeat(donor,repeats=total.shape[1]), #assign donor ID (SXX, X is a number)
                  "Library": np.repeat(library_number,repeats=total.shape[1]), # assing library number
                  "Sex": np.repeat(sex,repeats=total.shape[1]),
                  "Age": np.repeat(age,repeats=total.shape[1]),
                  "PMI_h": np.repeat(pmi_h,repeats=total.shape[1]),
                  "num_umi":np.repeat(num_umi,repeats=total.shape[1]),
                  "mean_reads_per_umi":np.repeat(mean_reads_per_umi,repeats=total.shape[1]),
                  "num_reads":np.repeat(num_reads,repeats=total.shape[1]),
                  "p_valid_barcodes":np.repeat(p_valid_barcodes,repeats=total.shape[1]),
                  "p_sequencing_saturation":np.repeat(p_sequencing_saturation,repeats=total.shape[1]),
                  "p_genome_not_gene":np.repeat(p_genome_not_gene,repeats=total.shape[1]),
                  "p_mapped_reads":np.repeat(p_mapped_reads,repeats=total.shape[1]),
                  "p_unmapped_reads":np.repeat(p_unmapped_reads,repeats=total.shape[1]),
                  "p_cell_del_filt":np.repeat(p_cell_del_filt,repeats=total.shape[1]),
                  "mean_counts_per_barcode":np.repeat(mean_counts_per_barcode,repeats=total.shape[1]),
                  "std_counts_per_barcode":np.repeat(std_counts_per_barcode,repeats=total.shape[1]),
                  "median_gpc":np.repeat(median_gpc,repeats=total.shape[1])}
    loompy.create(filename, total, row_attrs, col_attrs)
    
    print('restuctured loompy file ('+filename+') is saved!')
        #s.close()
    return

def get_rare_genes(D,donor,TH_min_number_of_cells_expr_a_gene):
    CPG = calc_cells_per_gene(D,donor)
    genes_to_remove = [CPG <= TH_min_number_of_cells_expr_a_gene]
    n_genes_discarded = np.sum(genes_to_remove[0]==True)
    return n_genes_discarded, genes_to_remove

def get_cells_with_many_genes(D,donor,TH_min_number_genes):
    GPC = calc_genes_per_cell(D,donor)
    cells_kept = [GPC >= TH_min_number_genes]
    n_cells_discarded = len(GPC)-np.sum(cells_kept)
    return n_cells_discarded, cells_kept
        
def get_cells_with_low_mito_fraction(D,donor,TH_max_fraction_mito_genes):
    counts_per_cell = calc_counts_per_cell(D,"",donor)
    gene_list = (D.ra.Gene).astype('str').tolist()
    id_mito_genes = [i for i,item in enumerate(gene_list) if item.startswith('MT-')]
    mito_counts_per_cell = calc_specific_counts_per_cell(id_mito_genes,D,donor)
    PMito = calc_percent_specific_counts(mito_counts_per_cell,counts_per_cell)
    cells_kept = [PMito/100 <= TH_max_fraction_mito_genes]
    n_cells_discarded = len(PMito)-np.sum(cells_kept)
    return n_cells_discarded, cells_kept
    
def get_cells_with_low_ribosomal_fraction(D,donor,TH_max_fraction_ribo_genes):
    counts_per_cell = calc_counts_per_cell(D,"",donor)
    gene_list = (D.ra.Gene).astype('str').tolist()
    #Rn45s or Rn4.5s
    id_ribo_genes = [i for i,item in enumerate(gene_list) if item.startswith('Rn4')]
    ribo_counts_per_cell = calc_specific_counts_per_cell(id_ribo_genes,D,donor)
    PRibo = calc_percent_specific_counts(ribo_counts_per_cell,counts_per_cell)
    cells_kept = [PRibo/100 <= TH_max_fraction_ribo_genes]
    n_cells_discarded = len(counts_per_cell)-np.sum(cells_kept)
    return n_cells_discarded, cells_kept

def get_cells_with_high_count_depth(D,donor,layer,TH_min_number_counts):
    counts_per_cell = calc_counts_per_cell(D,layer,donor)
    cells_kept = [counts_per_cell >= TH_min_number_counts]
    n_cells_discarded = len(counts_per_cell)-np.sum(cells_kept)
    return n_cells_discarded, cells_kept

def get_raw_data(library_nrs,samples_to_drop,opt_create_and_merge_raw_data_loom_files,path_project,path_code,path_data,path_summary_metrics):
    #merge raw data loom files if not already done
    if opt_create_and_merge_raw_data_loom_files ==True:
        path_patient_info_data = path_project+'/3_quality_control/data/'
        create_and_merge_loom_files(path_data,path_patient_info_data,path_summary_metrics,path_code,library_nrs,samples_to_drop)
    
    #get raw data
    os.chdir(path_data)
    D_raw = loompy.connect('Samples'+add_str+'.loom')
    os.chdir(path_code)
    return D_raw

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


def get_stats_per_sample(D,path_results,path_code,opt_stats_data,filtering_str,opt_save):

    if opt_stats_data == 'calculate':

        ### calculate some statistics per sample:
        # - CPG: cells/ counts per gene
        # - MCG: mean counts per gene
        # - min_CPG: minimum counts per gene
        # - GPC: genes per cell
        # - CPC: counts per cell
        # - MCPC: mitochondrial counts per cell
        # - PMito: percentage of counts belonging to mitochondrial mRNA per cell
        # - RCPC: ribosomal counts per cell
        # - PRibo: percentage of counts belonging to ribosomal mRNA per cell
        # - MeanCPC: mean counts per cell
        # - RLD: distribution of read length (1 value in distribution per gene)
        # - med_gpc: median genes per cell
        # - CPD: nuclei/ cells per donor
        # upcoming: 
        # - ProCodRPC: counts from protein coding genes (per cell)
        MCG, CPG, min_CPG, GPC, CPC, sCPC, usCPC, MCPC, PMito, RCPC, PRibo, mean_counts_per_cell, std_counts_per_cell, RLD, med_gpc, CPD = calc_stats_per_sample(D)

        ### Data frame with statistics:
        sample_IDs,n_cells,disease,sex = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
        
        df_stats = pd.DataFrame(data={'Donor ID': sample_IDs, 'Disease': disease, 'Sex': sex, 'Number of Cells': n_cells, 'Mean counts per barcode': mean_counts_per_cell, 'Std counts per barcode': std_counts_per_cell, 'Minimum number of counts per gene': min_CPG, 'Median genes per cell': med_gpc, 'Percentage mitochondrial counts per cell': PMito, 'Percentage ribosomal counts per cell': PRibo, 'Number of Cells':CPD })
    
        if opt_save == True:
            os.chdir(path_results)
            np.save('stats_per_sample'+filtering_str+'_MCG',MCG)
            np.save('stats_per_sample'+filtering_str+'_CPG',CPG)
            np.save('stats_per_sample'+filtering_str+'_min_CPG',min_CPG)
            np.save('stats_per_sample'+filtering_str+'_GPC',GPC)
            np.save('stats_per_sample'+filtering_str+'_CPC',CPC)
            np.save('stats_per_sample'+filtering_str+'_sCPC',sCPC)
            np.save('stats_per_sample'+filtering_str+'_usCPC',usCPC)
            np.save('stats_per_sample'+filtering_str+'_MCPC',MCPC)
            np.save('stats_per_sample'+filtering_str+'_PMito',PMito)
            np.save('stats_per_sample'+filtering_str+'_RCPC',RCPC)
            np.save('stats_per_sample'+filtering_str+'_PRibo',PRibo)
            np.save('stats_per_sample'+filtering_str+'_mean_counts_per_cell',mean_counts_per_cell)
            np.save('stats_per_sample'+filtering_str+'_std_counts_per_cell',std_counts_per_cell)
            np.save('stats_per_sample'+filtering_str+'_RLD',RLD)
            np.save('stats_per_sample'+filtering_str+'_median_genes_per_cell',med_gpc)
            #np.save('stats_per_sample'+filtering_str+'_cells_per_donor',CPD)
            #np.load(filename)
            df_stats.to_pickle('df_stats'+filtering_str)  # where to save it, usually as a .pkl
            # df = pd.read_pickle(file_name)
            os.chdir(path_code)
    elif opt_stats_data=='load':
        print(path_results)
        os.chdir(path_results)
        MCG = np.load('stats_per_sample'+filtering_str+'_MCG.npy',allow_pickle=True)
        CPG = np.load('stats_per_sample'+filtering_str+'_CPG.npy',allow_pickle=True)
        min_CPG = np.load('stats_per_sample'+filtering_str+'_min_CPG.npy',allow_pickle=True)
        GPC = np.load('stats_per_sample'+filtering_str+'_GPC.npy',allow_pickle=True)
        CPC = np.load('stats_per_sample'+filtering_str+'_CPC.npy',allow_pickle=True)
        sCPC = np.load('stats_per_sample'+filtering_str+'_sCPC.npy',allow_pickle=True)
        usCPC = np.load('stats_per_sample'+filtering_str+'_usCPC.npy',allow_pickle=True)
        MCPC = np.load('stats_per_sample'+filtering_str+'_MCPC.npy',allow_pickle=True)
        PMito = np.load('stats_per_sample'+filtering_str+'_PMito.npy',allow_pickle=True)
        RCPC = np.load('stats_per_sample'+filtering_str+'_MCPC.npy',allow_pickle=True)
        PRibo = np.load('stats_per_sample'+filtering_str+'_PRibo.npy',allow_pickle=True)
        mean_counts_per_cell = np.load('stats_per_sample'+filtering_str+'_mean_counts_per_cell.npy',allow_pickle=True)
        std_counts_per_cell = np.load('stats_per_sample'+filtering_str+'_std_counts_per_cell.npy',allow_pickle=True)
        RLD = np.load('stats_per_sample'+filtering_str+'_RLD.npy',allow_pickle=True)
        med_gpc = np.load('stats_per_sample'+filtering_str+'_median_genes_per_cell.npy',allow_pickle=True)
        #CPD = np.load('stats_per_sample'+filtering_str+'_cells_per_donor.npy',allow_pickle=True)
        df_stats = pd.read_pickle('df_stats'+filtering_str)
        os.chdir(path_code)
    else:
        print('The specified option for opt_stats_raw_data is not implemented!')
    return MCG, CPG, min_CPG, GPC, CPC, sCPC, usCPC, MCPC, PMito, RCPC, PRibo, mean_counts_per_cell, std_counts_per_cell, RLD, med_gpc, df_stats #CPD, df_stats

def get_filtering_str(opt_TH_filtered,opt_doublets_removed,opt_manual_doublet_TH_adjustment,opt_D_version,opt_pagoda_filtering):
    if opt_TH_filtered==True and opt_doublets_removed: 
        if opt_manual_doublet_TH_adjustment:
            filtering_str = '_TH_and_D_filtered_adj'+'_'+opt_D_version
        else:
            filtering_str = '_TH_and_D_filtered'
        if opt_pagoda_filtering:
            filtering_str = '_pagoda' + filtering_str
    elif opt_TH_filtered==True and opt_doublets_removed==False:
        filtering_str = '_TH_filtered'
        if opt_pagoda_filtering:
            filtering_str = '_pagoda' + filtering_str
    else:
        filtering_str='' 
    return filtering_str

def get_TH_filtered_data(sample_IDs, disease_status, path_code, path_data_raw, path_data_filt, path_results, opt_merge_raw_data_loom_files,opt_use_TH_filtered_data,opt_pagoda_filtering):
    if opt_use_TH_filtered_data == False:
        ## Quality control: Filtering
        ## create loom file with all filtered CTRL samples
        ## QC2: (1) remove genes that are rarely expressed: cells in which gene is expressed/ total #cells <= TH_min_fraction_of_cells_expr_a_gene
        ## QC2: (2) remove cells for which #genes <= TH_min_number_genes
        ## QC2: (4) remove cells with low count depth
        ## QC2: (3) remove cells for which number mitochondrial reads/ total number of reads >= TH_max_fraction_mito_genes
        df_TH_filtering_result = create_TH_filtered_loom_files(sample_IDs,disease_status,opt_pagoda_filtering,path_data_raw,path_data_filt,path_code,path_results)
  
    if opt_pagoda_filtering:
        add_str_pagoda = '_pagoda'
    else:
        add_str_pagoda = '_pagoda'
        
    os.chdir(path_data_filt)
    D_fil_CTRL = loompy.connect("Samples_CTRL"+add_str+add_str_pagoda+"_TH_filtered.loom")
    D_fil_SCZ = loompy.connect("Samples_SCZ"+add_str+add_str_pagoda+"_TH_filtered.loom")
    D_fil = loompy.connect("Samples"+add_str+add_str_pagoda+"_TH_filtered.loom")
    if opt_use_TH_filtered_data == True:
        os.chdir(path_results)
        df_TH_filtering_result = pd.read_pickle("df"+add_str_pagoda+"_TH_filtering_result"+add_str)
    os.chdir(path_code)
    return D_fil, D_fil_CTRL, D_fil_SCZ, df_TH_filtering_result

def get_TH_and_D_filtered_data(sample_IDs, disease_status, path_code, path_data_raw, path_data_filt, path_results, opt_pagoda_filtering, opt_merge_raw_data_loom_files, opt_use_TH_filtered_data, opt_use_TH_and_D_filtered_data, opt_manual_doublet_TH_adjustment, opt_plot_filtered_data, TH_min_cells_per_sample, opt_D_version):
        
    if opt_manual_doublet_TH_adjustment:
        add_str2 = "_adj"+'_'+opt_D_version
    else:
        add_str2 = ""
    
    if opt_pagoda_filtering:
        add_str_pagoda = '_pagoda'
    else:
        add_str_pagoda = '_pagoda'
        
    if opt_use_TH_and_D_filtered_data:
        os.chdir(path_data_filt)
        D_fil_CTRL = loompy.connect("Samples_CTRL"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        D_fil_SCZ = loompy.connect("Samples_SCZ"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        D_fil = loompy.connect("Samples"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        os.chdir(path_results)
        df_TH_and_doublet_filtering_result= pd.read_pickle("df"+add_str_pagoda+"_TH_and_D"+add_str2+"filtering_result"+add_str)
        os.chdir(path_data_filt)
        PD= np.load("PD"+add_str+add_str_pagoda+".npy",allow_pickle=True)
        #SCRUB= np.load("SCRUB"+add_str+".npy",allow_pickle=True)
        #DS= np.load("DS"+add_str+".npy",allow_pickle=True)
        #PDoublet= np.load("PDoublet"+add_str+".npy",allow_pickle=True)
        PD_updated= np.load("PD_updated"+add_str+add_str_pagoda+".npy",allow_pickle=True)
        os.chdir(path_code)
    else:
        D_TH_fil, _, _ , _ = get_TH_filtered_data(sample_IDs, disease_status, path_code, path_data_raw, path_data_filt, path_results, opt_merge_raw_data_loom_files,opt_use_TH_filtered_data,opt_pagoda_filtering)  
        TH_min_number_of_cells_expr_a_gene_all_samples, _, TH_min_number_counts, _, _, _, _, _,= get_TH_values()
        sample_IDs_TH_fil,_,disease_status_TH_fil,_ = get_sampleIDS_nCells_diseaseStatus_sexStatus(D_TH_fil)
        D_TH_fil.close()
        gc.collect()
        ## QC2: (3) remove doublets
        #doublet detection algorithm: scrublet, Wolock et al. 2019
        for d_id,donor in enumerate(sample_IDs_TH_fil):
            print(donor)
            scrub, doublet_scores, predicted_doublets = calc_doublet_prediction_scrublet(donor, add_str, add_str_pagoda, TH_min_number_of_cells_expr_a_gene_all_samples, TH_min_number_counts, path_data_filt, path_results, opt_plot_filtered_data)  
            if d_id==0:
                DS = [doublet_scores]
                PD = [predicted_doublets]
                SCRUB = [scrub]
                PDoublet = [np.sum(predicted_doublets)]
            else:
                DS.append(doublet_scores)
                PD.append(predicted_doublets)
                SCRUB.append(scrub)
                PDoublet.append(np.sum(predicted_doublets))
        os.chdir(path_data_filt)
        np.save('PD'+add_str+add_str_pagoda,PD)
        #np.save('SCRUB'+add_str,SCRUB)
        #np.save('DS'+add_str,DS)
        #np.save('PDoublet'+add_str,PDoublet)
        os.chdir(path_code)
        gc.collect()
        #adjust threshold manually
        if opt_manual_doublet_TH_adjustment:
            TH_doublet_score = get_TH_doublet_score(sample_IDs_TH_fil,opt_D_version)
            PD_updated = PD
            for d_id,donor in enumerate(sample_IDs_TH_fil):
                if TH_doublet_score[d_id]>0:
                    PD_updated[d_id] = SCRUB[d_id].call_doublets(threshold=TH_doublet_score[d_id])
            df_TH_and_doublet_filtering_result = create_TH_and_D_filtered_loom_files(sample_IDs_TH_fil,disease_status_TH_fil,PD_updated,opt_manual_doublet_TH_adjustment,opt_pagoda_filtering,add_str,TH_min_cells_per_sample,path_data_raw, path_data_filt,path_code,path_results,opt_D_version)
        else:
            PD_updated = []
            df_TH_and_doublet_filtering_result = create_TH_and_D_filtered_loom_files(sample_IDs_TH_fil,disease_status_TH_fil,PD,opt_manual_doublet_TH_adjustment,opt_pagoda_filtering,add_str,TH_min_cells_per_sample,path_data_raw, path_data_filt,path_code,path_results,opt_D_version)
        #load
        os.chdir(path_data_filt)
        np.save('PD_updated'+add_str+add_str_pagoda,PD_updated)
        D_fil_CTRL = loompy.connect("Samples_CTRL"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        D_fil_SCZ = loompy.connect("Samples_SCZ"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        D_fil = loompy.connect("Samples"+add_str+add_str_pagoda+"_TH_and_D"+add_str2+"_filtered.loom")
        os.chdir(path_code)
    return D_fil, D_fil_CTRL, D_fil_SCZ, df_TH_and_doublet_filtering_result, PD, PD_updated

def get_TH_doublet_score(sample_IDs,opt_D_version):
    TH_doublet_score = [0]*len(sample_IDs)
    if 'S1' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S1')[0][0]] = 0.2
    if 'S2' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S2')[0][0]] = 0.2   
    if 'S7' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S7')[0][0]] = 0.22   
    if 'S8' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S8')[0][0]] = 0.25  
    if 'S10' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S10')[0][0]] = 0.25
    if 'S11' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S11')[0][0]] = 0.25    
    if 'S12' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S12')[0][0]] = 0.25
    if 'S14' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S14')[0][0]] = 0.22
    if 'S16' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S16')[0][0]] = 0.21
    if 'S17' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S17')[0][0]] = 0.22
    if 'S18' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S18')[0][0]] = 0.22
    if 'S22' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S22')[0][0]] = 0.22    
    if 'S23' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S23')[0][0]] = 0.22   
    if 'S26' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S26')[0][0]] = 0.25
    if 'S27' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S27')[0][0]] = 0.25
    if 'S28' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S28')[0][0]] = 0.3
    if 'S31' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S31')[0][0]] = 0.25
    if 'S32' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S32')[0][0]] = 0.2  
    if 'S33' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S33')[0][0]] = 0.2 
    if 'S34' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S34')[0][0]] = 0.2
    if 'S35' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S35')[0][0]] = 0.2
    if 'S37' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S37')[0][0]] = 0.22
    if 'S38' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S38')[0][0]] = 0.2
    if 'S42' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S42')[0][0]] = 0.2
    if 'S49' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S49')[0][0]] = 0.25
    if 'S54' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S54')[0][0]] = 0.3
    if 'S64' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S64')[0][0]] = 0.3
    if 'S70' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S70')[0][0]] = 0.22
    if 'S73' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S73')[0][0]] = 0.2
    if 'S75' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S75')[0][0]] = 0.21
    if 'S80' in list(sample_IDs):
        TH_doublet_score[np.where(sample_IDs=='S80')[0][0]] = 0.21
                
    return TH_doublet_score

##plot:
def plot_CPC_vs_GPC(sample_IDs, disease, CPC, GPC, PM, PD, TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer_info, opt_colorCode,opt_save):
    plt = loadPltSettings(FS,10)
    color_def = cm.get_cmap('viridis', len(CPC))
    if opt_colorCode == 'mitoRPC_gradient':
        colmap = plt.cm.get_cmap('Blues')
        filename = layer_info+'_cd_vs_g_per_cell_colored_by_mito_counts_all_s'
    elif opt_colorCode == 'mitoRPC_cutoff':
        filename = layer_info+'_cd_vs_g_per_cell_colored_by_mito_counts_TH'
    elif opt_colorCode == 'doublet_prediction':
        filename = layer_info+'_cd_vs_g_per_cell_colored_by_doublets'
    elif opt_colorCode == 'group':
        filename = layer_info+'_cd_vs_g_per_cell_colored_by_group'
    elif opt_colorCode == 'sample':
        filename = layer_info+'_cd_vs_g_per_cell_colored_by_sample'
    else:
        filename = layer_info+'_cd_vs_g_per_cell'
    
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    
    #in case of scatter hist:
    if opt_colorCode == 'group' or opt_colorCode == 'sample':
        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.05
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.3]
        rect_histy = [left + width + spacing, bottom, 0.3, height]

    if opt_colorCode == 'sample':
        f6 = plt.figure(figsize=(8,8))  
        ax6 = f6.add_axes(rect_scatter)
        for d_id,donor in enumerate(sample_IDs):
            #sample 10% of cells for each donor:
            np.random.seed(1243431)
            idx = np.random.random_integers(0,len(CPC[d_id])-1,int(np.floor(len(CPC[d_id])*0.1)))
            ax6.scatter(CPC[d_id][idx],GPC[d_id][idx],s=5,facecolors=color_def(d_id),edgecolor='face',label=donor)
            ax6.set_xlabel('Counts per cell (log10 scale)')
            ax6.set_ylabel('Genes per cell (log10 scale)')
            ax6.set_yscale('log')
            ax6.set_xscale('log')
        ax6.legend(loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints=1)
        plt.show()
        if opt_save:
            f6.savefig(path_results+filename+'_colored_by_sample_all_samples_together.pdf', bbox_inches='tight')
        plt.close(f6)
        plt.clf()
    else:
        n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
        n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
        f6_CTRL, ax6_CTRL = plt.subplots(n_rows_CTRL,  n_cols, figsize=(35,30))
        f6_SCZ, ax6_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(35,30))
        if opt_colorCode == 'group':
            f6a = plt.figure(figsize=(8,8))  
            f6b = plt.figure(figsize=(8,8)) 
        else:
            f6, ax6 = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(35,12))
            f6.subplots_adjust(hspace = .15, wspace = 0.15)  
        f6_CTRL.subplots_adjust(hspace = .25, wspace = 0.25)  
        f6_SCZ.subplots_adjust(hspace = .25, wspace = 0.25)  
    
        c_id_CTRL=-1
        c_id_SCZ=-1
        r_id_CTRL=0
        r_id_SCZ=0
        for d_id,donor in enumerate(sample_IDs):
            if disease[d_id] =='CTRL':
                c_id_CTRL = c_id_CTRL + 1
                if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                    r_id_CTRL = r_id_CTRL + 1
                    c_id_CTRL = 0
                v_max=np.max(PM[d_id])
                if opt_colorCode == 'mitoRPC_gradient':
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id],GPC[d_id],c=PM[d_id],vmin=0,vmax=v_max,s=3,cmap=colmap,edgecolors='face',alpha=0.7) 
                elif opt_colorCode == 'mitoRPC_cutoff':
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id][PM[d_id]<=TH_max_fraction_mito_counts*100],GPC[d_id][PM[d_id]<=TH_max_fraction_mito_counts*100],color='grey',s=3,label="PM<=5%" if ((r_id_CTRL == (n_rows_CTRL-1)) & (c_id_CTRL == (n_cols-1))) else '') 
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id][PM[d_id]>TH_max_fraction_mito_counts*100],GPC[d_id][PM[d_id]>TH_max_fraction_mito_counts*100],color='blue',s=3,label='PM>5%' if ((r_id_CTRL == (n_rows_CTRL-1)) & (c_id_CTRL == (n_cols-1))) else '') 
                elif opt_colorCode == 'doublet_prediction':
                    #To Do: needs to be fixed. np.shape(GPC_f[0]) is not equal to np.shape(PD[0]) or np.shape(PD_updated[0])
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id],GPC[d_id],color='grey',s=3,label='predicted single nuclei' if ((r_id_CTRL == (n_rows_CTRL-1)) & (c_id_CTRL == (n_cols-1))) else '') 
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id][PD[d_id]],GPC[d_id][PD[d_id]],color='red',s=3,label='predicted doublet' if ((r_id_CTRL == (n_rows_CTRL-1)) & (c_id_CTRL == (n_cols-1))) else '')
                elif opt_colorCode == 'group':
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id][np.logical_or(CPC[d_id]>=TH_min_number_counts,GPC[d_id]>=TH_min_number_genes)],GPC[d_id][np.logical_or(CPC[d_id]>=TH_min_number_counts, GPC[d_id]>=TH_min_number_genes)],s=2,color = 'steelblue') 
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id][np.logical_or(CPC[d_id]<TH_min_number_counts,GPC[d_id]<TH_min_number_genes)],GPC[d_id][np.logical_or(CPC[d_id]<TH_min_number_counts,GPC[d_id]<TH_min_number_genes)],s=2,color = 'lightgrey') 
                else:
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].scatter(CPC[d_id],GPC[d_id],s=2) 
                ax6_CTRL[r_id_CTRL,c_id_CTRL].title.set_text(donor+' ('+disease[d_id][0]+')')
                ax6_CTRL[r_id_CTRL,c_id_CTRL].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                xlim = ax6_CTRL[r_id_CTRL,c_id_CTRL].get_xlim()
                ylim = ax6_CTRL[r_id_CTRL,c_id_CTRL].get_ylim()
                ax6_CTRL[r_id_CTRL,c_id_CTRL].plot(xlim,[TH_min_number_genes,TH_min_number_genes],'-r',linewidth=1)
                ax6_CTRL[r_id_CTRL,c_id_CTRL].plot([TH_min_number_counts,TH_min_number_counts],ylim,'-r',linewidth=1)
                #ax[r_id,c_id].plot([TH_min_counts[d_id],TH_min_counts[d_id]], [0, np.max(genes_per_cell)],'--r',linewidth=0.5)
                if opt_colorCode == 'mitoRPC_gradient':
                    norm = colors.Normalize(vmin=0, vmax=v_max)
                    sm =  cm.ScalarMappable(norm=norm, cmap=colmap)
                    sm.set_array([])
                    f6_CTRL.colorbar(sm, ax=ax6_CTRL[r_id_CTRL,c_id_CTRL])
                    #cbar = f6_CTRL.colorbar(sm, ax=ax6_CTRL[r_id_CTRL,c_id_CTRL])
                    #cbar.ax.set_title("PM")
                elif opt_colorCode == 'mitoRPC_cutoff' or opt_colorCode == 'doublet_prediction':
                    ax6_CTRL[r_id_CTRL,c_id_CTRL].legend()
            else:
                c_id_SCZ = c_id_SCZ + 1
                if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                    r_id_SCZ = r_id_SCZ + 1
                    c_id_SCZ = 0
                v_max=np.max(PM[d_id])
                if opt_colorCode == 'mitoRPC_gradient':
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id],GPC[d_id],c=PM[d_id],vmin=0,vmax=v_max,s=3,cmap=colmap,edgecolors='face',alpha=0.7) 
                elif opt_colorCode == 'mitoRPC_cutoff':
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id][PM[d_id]<=TH_max_fraction_mito_counts*100],GPC[d_id][PM[d_id]<=TH_max_fraction_mito_counts*100],color='grey',s=3,label='PM<=5%' if ((r_id_SCZ == n_rows_SCZ-1) and (c_id_SCZ == n_cols-1)) else '') 
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id][PM[d_id]>TH_max_fraction_mito_counts*100],GPC[d_id][PM[d_id]>TH_max_fraction_mito_counts*100],color='blue',s=3,label='PM>5%' if ((r_id_SCZ == n_rows_SCZ-1) and (c_id_SCZ == n_cols-1)) else '') 
                elif opt_colorCode == 'doublet_prediction':
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id],GPC[d_id],color='grey',s=3,label='predicted single nuclei' if ((r_id_SCZ == n_rows_SCZ-1) and (c_id_SCZ == n_cols-1)) else '') 
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id][PD[d_id]],GPC[d_id][PD[d_id]],color='red',s=3,label='predicted doublet' if ((r_id_SCZ == n_rows_SCZ-1) and (c_id_SCZ == n_cols-1)) else '')
                elif opt_colorCode == 'group':
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id][np.logical_or(CPC[d_id]>=TH_min_number_counts, GPC[d_id]>=TH_min_number_genes)],GPC[d_id][np.logical_or(CPC[d_id]>=TH_min_number_counts, GPC[d_id]>=TH_min_number_genes)],s=2,color='darkorange') 
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id][np.logical_or(CPC[d_id]<TH_min_number_counts, GPC[d_id]<TH_min_number_genes)],GPC[d_id][np.logical_or(CPC[d_id]<TH_min_number_counts, GPC[d_id]<TH_min_number_genes)],s=2,color='lightgrey') 
                else:
                    ax6_SCZ[r_id_SCZ,c_id_SCZ].scatter(CPC[d_id],GPC[d_id],s=2) 
                ax6_SCZ[r_id_SCZ,c_id_SCZ].title.set_text(donor+' ('+disease[d_id][0]+')')
                ax6_SCZ[r_id_SCZ,c_id_SCZ].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                xlim = ax6_SCZ[r_id_SCZ,c_id_SCZ].get_xlim()
                ylim = ax6_SCZ[r_id_SCZ,c_id_SCZ].get_ylim()
                ax6_SCZ[r_id_SCZ,c_id_SCZ].plot(xlim,[TH_min_number_genes,TH_min_number_genes],'-r',linewidth=1)
                ax6_SCZ[r_id_SCZ,c_id_SCZ].plot([TH_min_number_counts,TH_min_number_counts],ylim,'-r',linewidth=1)
                #ax[r_id,c_id].plot([TH_min_counts[d_id],TH_min_counts[d_id]], [0, np.max(genes_per_cell)],'--r',linewidth=0.5)
                if opt_colorCode == 'mitoRPC_gradient':
                    norm = colors.Normalize(vmin=0, vmax=v_max)
                    sm =  cm.ScalarMappable(norm=norm, cmap=colmap)
                    sm.set_array([])
                    f6_SCZ.colorbar(sm, ax=ax6_SCZ[r_id_SCZ,c_id_SCZ])
                    #cbar = f6_SCZ.colorbar(sm, ax=ax6_SCZ[r_id_SCZ,c_id_SCZ])
                    #cbar.ax.set_title("PM")
                elif opt_colorCode == 'mitoRPC_cutoff' or opt_colorCode == 'doublet_prediction':
                        ax6_SCZ[r_id_SCZ,c_id_SCZ].legend()
        colors_groups = ['steelblue','darkorange']
        for g_id in range(0,len(G)):
            d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
            if opt_colorCode == 'mitoRPC_gradient':
                ax6[g_id].scatter(np.hstack(CPC[d_IDs]),np.hstack(GPC[d_IDs]),c=np.hstack(PM[d_IDs]),vmin=0,vmax=v_max,s=3,cmap=colmap,edgecolors='face',alpha=0.7) 
                norm = colors.Normalize(vmin=0, vmax=v_max)
                sm =  cm.ScalarMappable(norm=norm, cmap=colmap)
                sm.set_array([])
                f6.colorbar(sm, ax=ax6[g_id])
                #cbar = f6.colorbar(sm, ax=ax6[g_id])
                #cbar.ax.set_title("PM")
            elif opt_colorCode == 'mitoRPC_cutoff':
                ax6[g_id].scatter(np.hstack(CPC[d_IDs])[np.hstack(PM[d_IDs])<=TH_max_fraction_mito_counts*100],np.hstack(GPC[d_IDs])[np.hstack(PM[d_IDs])<=TH_max_fraction_mito_counts*100],color='grey',s=3,label='PM<=5%' if (g_id==1) else '') 
                ax6[g_id].scatter(np.hstack(CPC[d_IDs])[np.hstack(PM[d_IDs])>TH_max_fraction_mito_counts*100],np.hstack(GPC[d_IDs])[np.hstack(PM[d_IDs])>TH_max_fraction_mito_counts*100],color='blue',s=3,label='PM>5%' if (g_id==1) else '') 
                ax6_SCZ[r_id_SCZ,c_id_SCZ].legend()
            elif opt_colorCode == 'doublet_prediction':
                ax6[g_id].scatter(np.hstack(CPC[d_IDs]),np.hstack(GPC[d_IDs]),color='grey',s=3,label='predicted single nuclei' if (g_id==1) else '') 
                ax6[g_id].scatter(np.hstack(CPC[d_IDs])[np.hstack(PD[d_IDs])],np.hstack(GPC[d_IDs])[np.hstack(PD[d_IDs])],color='red',s=3,label='predicted doublet' if (g_id==1) else '')
            if opt_colorCode == 'group':
                if g_id==0:
                    ax6a = f6a.add_axes(rect_scatter)
                    #add GPC histogram to y axis
                    ax6a_histx = f6a.add_axes(rect_histx,sharex=ax6a)
                    if layer_info=="unspliced":
                        ax6a_histy = f6a.add_axes(rect_histy,sharey=ax6a)
                    # no labels
                    x=np.hstack(CPC[d_IDs])
                    y=np.hstack(GPC[d_IDs])
                    ax6a_histx.tick_params(axis="x", labelbottom=False)
                    if layer_info=="unspliced":
                        ax6a_histy.tick_params(axis="y", labelleft=False)
                        pos=0.3
                    else:
                        pos=0.4
                    ax6a.scatter(x[np.logical_or(x>=TH_min_number_counts, y>=TH_min_number_genes)],y[np.logical_or(x>=TH_min_number_counts, y>=TH_min_number_genes)],s=1,color = colors_groups[g_id])
                    ax6a.scatter(x[np.logical_or(x<TH_min_number_counts, y<TH_min_number_genes)],y[np.logical_or(x<TH_min_number_counts, y<TH_min_number_genes)],s=1,color = 'lightgrey')
                    xlim = ax6a.get_xlim()
                    ylim = ax6a.get_ylim()
                    ax6a.plot(xlim,[TH_min_number_genes,TH_min_number_genes],'-k',linewidth=2)
                    ax6a.plot([TH_min_number_counts,TH_min_number_counts],ylim,'-k',linewidth=2)
                    ax6a_histx.hist(x,color='darkgrey', bins=50)
                    if layer_info=="unspliced":
                        ax6a_histy.hist(y, color='darkgrey', orientation='horizontal', bins=50)
                    f6a.text(pos, 0.001, layer_info+' counts per cell', ha='center',fontsize = FS+10)
                    f6a.text(0.01, 0.6, 'Genes per cell', va='center', rotation='vertical',fontsize = FS+10)
                    plt.show()
                    if opt_save:
                        f6a.savefig(path_results+filename+'_'+G[g_id]+'_all_samples_together.pdf', bbox_inches='tight')
                    plt.close(f6a)
                    plt.clf()
                else:
                    ax6b = f6b.add_axes(rect_scatter)
                    #add GPC histogram to y axis
                    ax6b_histx = f6b.add_axes(rect_histx,sharex=ax6a)
                    if layer_info=="unspliced":
                        ax6b_histy = f6b.add_axes(rect_histy,sharey=ax6a)
                    # no labels
                    x=np.hstack(CPC[d_IDs])
                    y=np.hstack(GPC[d_IDs])
                    ax6b_histx.tick_params(axis="x", labelbottom=False)
                    if layer_info=="unspliced":
                        ax6b_histy.tick_params(axis="y", labelleft=False)
                        pos=0.4
                    else:
                        pos=0.4
                    ax6b.scatter(x[np.logical_or(x>=TH_min_number_counts, y>=TH_min_number_genes)],y[np.logical_or(x>=TH_min_number_counts, y>=TH_min_number_genes)],s=1,color = colors_groups[g_id])
                    ax6b.scatter(x[np.logical_or(x<TH_min_number_counts, y<TH_min_number_genes)],y[np.logical_or(x<TH_min_number_counts, y<TH_min_number_genes)],s=1,color = 'lightgrey')
                    xlim = ax6b.get_xlim()
                    ylim = ax6b.get_ylim()
                    ax6b.plot(xlim,[TH_min_number_genes,TH_min_number_genes],'-k',linewidth=2)
                    ax6b.plot([TH_min_number_counts,TH_min_number_counts],ylim,'-k',linewidth=2)
                    ax6b_histx.hist(x,color='darkgrey', bins=50)
                    if layer_info=="unspliced":
                        ax6b_histy.hist(y, color='darkgrey', orientation='horizontal', bins=50)
                    f6b.text(pos, 0.001, layer_info+' counts per cell', ha='center',fontsize = FS+10)
                    f6b.text(0.01, 0.6, 'Genes per cell', va='center', rotation='vertical',fontsize = FS+10)
                    plt.show()
                    if opt_save:
                        f6b.savefig(path_results+filename+'_'+G[g_id]+'_all_samples_together.pdf', bbox_inches='tight')
                    plt.close(f6b)
                    plt.clf()
            else:
                ax6[g_id].scatter(np.hstack(CPC[d_IDs]),np.hstack(GPC[d_IDs]),s=2) 
            if opt_colorCode != 'group':
                ax6[g_id].title.set_text('All samples '+G[g_id])
                ax6[g_id].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
                xlim = ax6[g_id].get_xlim()
                ylim = ax6[g_id].get_ylim()
                ax6[g_id].plot(xlim,[TH_min_number_genes,TH_min_number_genes],'-k',linewidth=1)
                ax6[g_id].plot([TH_min_number_counts,TH_min_number_counts],ylim,'-k',linewidth=1)
            
        if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
            for r_id_del in range(r_id_CTRL,n_rows_CTRL):
                for c_id_del in range(c_id_CTRL+1,n_cols):
                    f6_CTRL.delaxes(ax6_CTRL[r_id_del,c_id_del])
        if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
            for r_id_del in range(r_id_SCZ,n_rows_SCZ):
                for c_id_del in range(c_id_SCZ+1,n_cols):
                    f6_SCZ.delaxes(ax6_SCZ[r_id_del,c_id_del])
        f6_CTRL.text(0.5, 0.04, layer_info+' count depth ('+layer_info+' counts per cellID)', ha='center',fontsize = FS+10)
        f6_SCZ.text(0.5, 0.04, layer_info+' count depth ('+layer_info+' counts per cellID)', ha='center',fontsize = FS+10)
        if opt_colorCode != 'group':
            f6.text(0.5, 0.01, layer_info+' count depth ('+layer_info+' counts per cellID)', ha='center',fontsize = FS+10)
            f6.text(0.05, 0.5, 'Genes per cell ID', va='center', rotation='vertical',fontsize = FS+10)
            f6.suptitle(layer_info+' count depth vs. genes per barcode',fontsize = FS+10)
        f6_CTRL.text(0.05, 0.5, 'Genes per cell ID', va='center', rotation='vertical',fontsize = FS+10)
        f6_SCZ.text(0.05, 0.5, 'Genes per cell ID', va='center', rotation='vertical',fontsize = FS+10)
        f6_CTRL.suptitle(layer_info+' count depth vs. genes per barcode',fontsize = FS+10)
        f6_SCZ.suptitle(layer_info+' count depth vs. genes per barcode',fontsize = FS+10)
        if opt_save:
            f6_CTRL.savefig(path_results+filename+'_CTRL.pdf', bbox_inches='tight')
            f6_SCZ.savefig(path_results+filename+'_SCZ.pdf', bbox_inches='tight')
            if opt_colorCode != 'group':
                f6.savefig(path_results+filename+'_all.pdf', bbox_inches='tight')
            #f6.savefig(path_results+filename+'_all_samples.svg', bbox_inches='tight')
        if opt_colorCode != 'group':
            plt.close(f6)
        plt.close(f6_SCZ)
        plt.close(f6_CTRL)
        plt.clf()

def plot_RLD(sample_IDs, disease, RLD, path_results, opt_save):
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
    n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
    plt = loadPltSettings(FS,10)
    f7_CTRL, ax7_CTRL = plt.subplots(n_rows_CTRL, n_cols, figsize=(35,30))
    f7_SCZ, ax7_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(35,30))
    f7, ax7 = plt.subplots(1, 2, sharey=True, figsize=(35,12))
    f7_CTRL.subplots_adjust(hspace = .35, wspace = 0.15)  
    f7_SCZ.subplots_adjust(hspace = .35, wspace = 0.15)  
    f7.subplots_adjust(hspace = .35, wspace = 0.15)  
    c_id_CTRL=-1
    c_id_SCZ=-1
    r_id_CTRL=0
    r_id_SCZ=0
    for d_id,donor in enumerate(sample_IDs):
        if disease[d_id] =='CTRL':
            c_id_CTRL = c_id_CTRL + 1
            if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                r_id_CTRL = r_id_CTRL + 1
                c_id_CTRL = 0
            ax7_CTRL[r_id_CTRL,c_id_CTRL].hist(RLD[d_id],bins=100,density=True,color='steelblue')
        else:
            c_id_SCZ = c_id_SCZ + 1
            if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                r_id_SCZ = r_id_SCZ + 1
                c_id_SCZ = 0
            ax7_SCZ[r_id_SCZ,c_id_SCZ].hist(RLD[d_id],bins=100,density=True,color='steelblue')
            
    for g_id in range(0,len(G)):
        d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
        ax7[g_id].hist(np.hstack(RLD[d_IDs]),bins=100,density=True,color='steelblue')
    f7_CTRL.text(0.5, 0.04, 'Count depth (counts per cellID)', ha='center')
    f7_SCZ.text(0.5, 0.04, 'Count depth (counts per cellID)', ha='center')
    f7.text(0.5, 0.04, 'Count depth (counts per cellID)', ha='center')
    f7_CTRL.text(0.01, 0.5, 'Absolute frequency', va='center', rotation='vertical')
    f7_SCZ.text(0.01, 0.5, 'Absolute frequency', va='center', rotation='vertical')
    f7.text(0.01, 0.5, 'Absolute frequency', va='center', rotation='vertical')
    f7_CTRL.suptitle('Read length in bp',fontsize = FS+5)
    f7_SCZ.suptitle('Read length in bp',fontsize = FS+5)
    f7.suptitle('Read length in bp',fontsize = FS+5)      
    # if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
    #     for r_id_del in range(r_id_CTRL,n_rows_CTRL):
    #         for c_id_del in range(c_id_CTRL,n_cols):
    #             f7_CTRL.delaxes(ax7_CTRL[r_id_del,c_id_del])
    # if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
    #     for r_id_del in range(r_id_SCZ,n_rows_SCZ):
    #         for c_id_del in range(c_id_SCZ+1,n_cols):
    #             f7_SCZ.delaxes(ax7_SCZ[r_id_del,c_id_del])
    if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
        for r_id_del in range(r_id_CTRL,n_rows_CTRL):
            for c_id_del in range(c_id_CTRL+1,n_cols):
                f7_CTRL.delaxes(ax7_CTRL[r_id_del,c_id_del])
    if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
        for r_id_del in range(r_id_SCZ,n_rows_SCZ):
            for c_id_del in range(c_id_SCZ+1,n_cols):
                f7_SCZ.delaxes(ax7_SCZ[r_id_del,c_id_del])
    if opt_save: 
        f7_CTRL.savefig(path_results+'rl_distr_CTRL.pdf', bbox_inches='tight')
        f7_SCZ.savefig(path_results+'rl_distr_SCZ.pdf', bbox_inches='tight')
        f7.savefig(path_results+'rl_distr_all.pdf', bbox_inches='tight')
    plt.close(f7)
    plt.close(f7_SCZ)
    plt.close(f7_CTRL)
    plt.clf()

def plot_MCG(sample_IDs, disease, MCG, path_results, opt_save):
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
    n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
    plt = loadPltSettings(FS,10)
    f2_CTRL, ax2_CTRL = plt.subplots(n_rows_CTRL, n_cols, figsize=(37,30))
    f2_SCZ, ax2_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(37,30))
    f2, ax2 = plt.subplots(1, 2, sharey=True, figsize=(35,6))
    f2_CTRL.subplots_adjust(hspace = .35, wspace = 0.15)  
    f2_SCZ.subplots_adjust(hspace = .35, wspace = 0.15)  
    f2.subplots_adjust(hspace = .35, wspace = 0.15)  
    c_id_CTRL=-1
    c_id_SCZ=-1
    r_id_CTRL=0
    r_id_SCZ=0
    for d_id,donor in enumerate(sample_IDs):
        if disease[d_id] =='CTRL':
            c_id_CTRL = c_id_CTRL + 1
            if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                r_id_CTRL = r_id_CTRL + 1
                c_id_CTRL = 0
            ax2_CTRL[r_id_CTRL,c_id_CTRL].hist(MCG[d_id],bins=100,density=True,color='steelblue')
            ax2_CTRL[r_id_CTRL,c_id_CTRL].title.set_text(donor+' ('+disease[d_id]+')')
        else:
            c_id_SCZ = c_id_SCZ + 1
            if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                r_id_SCZ = r_id_SCZ + 1
                c_id_SCZ = 0
            ax2_SCZ[r_id_SCZ,c_id_SCZ].hist(MCG[d_id],bins=100,density=True,color='steelblue')
            ax2_SCZ[r_id_SCZ,c_id_SCZ].title.set_text(donor+' ('+disease[d_id]+')')
    for g_id in range(0,len(G)): 
        d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
        ax2[g_id].hist(MCG[d_IDs].flatten(),bins=100,density=True,color='steelblue')
        ax2[g_id].title.set_text('All donors belonging to '+G[g_id]+' group')
    f2_CTRL.text(0.5, 0.04, 'Average counts per gene', ha='center')
    f2_SCZ.text(0.5, 0.04, 'Average counts per gene', ha='center')
    f2.text(0.5, 0.02, 'Average counts per gene', ha='center')
    f2_CTRL.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    f2_SCZ.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    f2.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    # if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
    #     for r_id_del in range(r_id_CTRL,n_rows_CTRL):
    #         for c_id_del in range(c_id_CTRL,n_cols):
    #             f2_CTRL.delaxes(ax2_CTRL[r_id_del,c_id_del])
    # if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
    #     for r_id_del in range(r_id_SCZ,n_rows_SCZ):
    #         for c_id_del in range(c_id_SCZ+1,n_cols):
    #             f2_SCZ.delaxes(ax2_SCZ[r_id_del,c_id_del])
                
    if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
        for r_id_del in range(r_id_CTRL,n_rows_CTRL):
            for c_id_del in range(c_id_CTRL+1,n_cols):
                f2_CTRL.delaxes(ax2_CTRL[r_id_del,c_id_del])
    if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
        for r_id_del in range(r_id_SCZ,n_rows_SCZ):
            for c_id_del in range(c_id_SCZ+1,n_cols):
                f2_SCZ.delaxes(ax2_SCZ[r_id_del,c_id_del])
    f2_CTRL.suptitle('Distribution of mean counts per gene',fontsize = FS+5)
    f2_SCZ.suptitle('Distribution of mean counts per gene',fontsize = FS+5)
    f2.suptitle('Distribution of mean counts per gene',fontsize = FS+5)
    if opt_save:
        f2_CTRL.savefig(path_results+'mean_counts_per_gene_all_CTRL_samples.pdf', bbox_inches='tight')
        f2_SCZ.savefig(path_results+'mean_counts_per_gene_all_SCZ_samples.pdf', bbox_inches='tight')
        f2.savefig(path_results+'mean_counts_per_gene_all_samples.pdf', bbox_inches='tight')
    plt.close(f2_CTRL)
    plt.close(f2_SCZ)
    plt.close(f2)
    plt.clf()
    
def plot_GPC(sample_IDs, disease, sex, GPC, TH_min_number_genes, path_results, opt_save):
    color_spec=['steelblue','darkorange']
    color_spec2=['lightseagreen','darkorchid']
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
    n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
    plt = loadPltSettings(FS,10)
    f3_CTRL, ax3_CTRL = plt.subplots(n_rows_CTRL, n_cols, figsize=(37,30))
    f3_SCZ, ax3_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(37,30))
    f3, ax3 = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(35,12))
    f3b, ax3b = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(35,12))
    f3_CTRL.subplots_adjust(hspace = .5, wspace = 0.25)  
    f3_SCZ.subplots_adjust(hspace = .5, wspace = 0.25)  
    f3.subplots_adjust(hspace = .35, wspace = 0.25)  
    f3b.subplots_adjust(hspace = .35, wspace = 0.25)  
    c_id_CTRL=-1
    c_id_SCZ=-1
    r_id_CTRL=0
    r_id_SCZ=0
    #min_genes_per_cell = np.empty(len(sample_IDs), dtype=float)
    for d_id,donor in enumerate(sample_IDs):
        #min_genes_per_cell[d_id] = np.min(GPC[d_id])
        if disease[d_id] =='CTRL':
            c_id_CTRL = c_id_CTRL + 1
            if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                r_id_CTRL = r_id_CTRL + 1
                c_id_CTRL = 0
            ax3_CTRL[r_id_CTRL,c_id_CTRL].hist(GPC[d_id],bins=100,density=True,color='steelblue')
            ylim = ax3_CTRL[r_id_CTRL,c_id_CTRL].get_ylim()
            ax3_CTRL[r_id_CTRL,c_id_CTRL].plot([TH_min_number_genes, TH_min_number_genes], ylim,'-b')
            ax3_CTRL[r_id_CTRL,c_id_CTRL].title.set_text(donor+' ('+disease[d_id]+')')
            ax3_CTRL[r_id_CTRL,c_id_CTRL].ticklabel_format(style='scientific')
        else:
            c_id_SCZ = c_id_SCZ + 1
            if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                r_id_SCZ = r_id_SCZ + 1
                c_id_SCZ = 0
            ax3_SCZ[r_id_SCZ,c_id_SCZ].hist(GPC[d_id],bins=100,density=True,color='darkorange')
            ylim = ax3_SCZ[r_id_SCZ,c_id_SCZ].get_ylim()
            ax3_SCZ[r_id_SCZ,c_id_SCZ].plot([TH_min_number_genes, TH_min_number_genes], ylim,'-b')
            ax3_SCZ[r_id_SCZ,c_id_SCZ].title.set_text(donor+' ('+disease[d_id]+')')
            ax3_SCZ[r_id_SCZ,c_id_SCZ].ticklabel_format(style='scientific')
    for g_id in range(0,len(G)): 
        d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
        ax3[g_id].hist(np.hstack(GPC[d_IDs]),bins=100,density=True,color=color_spec[g_id])
        ax3[g_id].title.set_text('All donors belonging to '+G[g_id]+' group')
        ylim = ax3[g_id].get_ylim()
        ax3[g_id].plot([TH_min_number_genes, TH_min_number_genes], ylim,'-b')
    for g_id in range(0,len(G2)): 
        d_IDs = [i for i in range(0,len(sex)) if sex[i]==G2[g_id]]
        ax3b[g_id].hist(np.hstack(GPC[d_IDs]),bins=100,density=True,color=color_spec2[g_id])
        ax3b[g_id].title.set_text('All donors belonging to '+G2[g_id]+' group')
        ylim = ax3b[g_id].get_ylim()
        ax3b[g_id].plot([TH_min_number_genes, TH_min_number_genes], ylim,'-b')
    f3_CTRL.text(0.5, 0.01, 'Genes per cell ID', ha='center')
    f3_SCZ.text(0.5, 0.01, 'Genes per cell ID', ha='center')
    f3.text(0.5, 0.01, 'Genes per cell ID', ha='center')
    f3b.text(0.5, 0.01, 'Genes per cell ID', ha='center')
    f3_CTRL.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    f3_SCZ.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    f3.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    f3b.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
    # if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
    #     for r_id_del in range(r_id_CTRL,n_rows_CTRL):
    #         for c_id_del in range(c_id_CTRL,n_cols):
    #             f3_CTRL.delaxes(ax3_CTRL[r_id_del,c_id_del])
    # if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
    #     for r_id_del in range(r_id_SCZ,n_rows_SCZ):
    #         for c_id_del in range(c_id_SCZ+1,n_cols):
    #             f3_SCZ.delaxes(ax3_SCZ[r_id_del,c_id_del])
    if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
        for r_id_del in range(r_id_CTRL,n_rows_CTRL):
            for c_id_del in range(c_id_CTRL+1,n_cols):
                f3_CTRL.delaxes(ax3_CTRL[r_id_del,c_id_del])
    if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
        for r_id_del in range(r_id_SCZ,n_rows_SCZ):
            for c_id_del in range(c_id_SCZ+1,n_cols):
                f3_SCZ.delaxes(ax3_SCZ[r_id_del,c_id_del])
    f3_CTRL.suptitle('Distribution of genes per cellID (barcode)',fontsize = FS+5)
    f3_SCZ.suptitle('Distribution of genes per cellID (barcode)',fontsize = FS+5)
    f3.suptitle('Distribution of genes per cellID (barcode)',fontsize = FS+5)
    f3b.suptitle('Distribution of genes per cellID (barcode)',fontsize = FS+5)
    if opt_save:
        f3_CTRL.savefig(path_results+'genes_per_cell_distribution_all_CTRL_samples.pdf', bbox_inches='tight')
        f3_SCZ.savefig(path_results+'genes_per_cell_distribution_all_SCZ_samples.pdf', bbox_inches='tight')
        f3.savefig(path_results+'genes_per_cell_distribution_all_samples_grouped_by_disease_status.pdf', bbox_inches='tight')
        f3b.savefig(path_results+'genes_per_cell_distribution_all_samples_grouped_by_sex_status.pdf', bbox_inches='tight')
        #f3_CTRL.savefig(path_results+'genes_per_cell_distribution_all_CTRL_samples.svg', bbox_inches='tight')
        #f3_SCZ.savefig(path_results+'genes_per_cell_distribution_all_SCZ_samples.svg', bbox_inches='tight')
        #f3.savefig(path_results+'genes_per_cell_distribution_all_samples.svg', bbox_inches='tight')
        #f3b.savefig(path_results+'genes_per_cell_distribution_all_samples.svg', bbox_inches='tight')
    plt.close(f3)
    plt.close(f3b)
    plt.close(f3_SCZ)
    plt.close(f3_CTRL)
    plt.clf()

def plot_CPC(sample_IDs, disease, CPC, path_results, TH_min_counts, layer_info, opt_style, opt_save):
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
    n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
    plt = loadPltSettings(FS,10)
    f4_CTRL, ax4_CTRL = plt.subplots(n_rows_CTRL, n_cols,figsize=(35,30))
    f4_SCZ, ax4_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(35,30))
    if opt_style == 'dots':
        f4, ax4 = plt.subplots(1, 2, sharey=True, sharex=True, figsize=(35,25))
    else:
        f4, ax4 = plt.subplots(1, 2, sharey=True, figsize=(35,25))
    f4_CTRL.subplots_adjust(hspace = .5, wspace = 0.15)
    f4_SCZ.subplots_adjust(hspace = .5, wspace = 0.15)  
    f4.subplots_adjust(hspace = .35, wspace = 0.15)  
    c_id_CTRL=-1
    c_id_SCZ=-1
    r_id_CTRL=0
    r_id_SCZ=0
    for d_id,donor in enumerate(sample_IDs):
        #min_genes_per_cell[d_id] = np.min(GPC[d_id])
        if disease[d_id] =='CTRL':
            c_id_CTRL = c_id_CTRL + 1
            if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                r_id_CTRL = r_id_CTRL + 1
                c_id_CTRL = 0
            if opt_style == 'histogram':
                ax4_CTRL[r_id_CTRL,c_id_CTRL].hist(CPC[d_id],bins=100,density=True,color='steelblue')
                ylim = ax4_CTRL[r_id_CTRL,c_id_CTRL].get_ylim()
                ax4_CTRL[r_id_CTRL,c_id_CTRL].plot([TH_min_counts, TH_min_counts], ylim,'-k')
            elif opt_style == 'dots':
                ax4_CTRL[r_id_CTRL,c_id_CTRL].scatter(range(1,len(CPC[d_id])+1),np.sort(CPC[d_id]),color = 'steelblue',s=3)
                xlim = ax4_CTRL[r_id_CTRL,c_id_CTRL].get_xlim()
                ax4_CTRL[r_id_CTRL,c_id_CTRL].plot(xlim,[TH_min_counts, TH_min_counts], '-k')
            ax4_CTRL[r_id_CTRL,c_id_CTRL].title.set_text(donor+' ('+disease[d_id]+')')
        else:
            c_id_SCZ = c_id_SCZ + 1
            if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                r_id_SCZ = r_id_SCZ + 1
                c_id_SCZ = 0
            if opt_style == 'histogram':
                ax4_SCZ[r_id_SCZ,c_id_SCZ].hist(CPC[d_id],bins=100,density=True,color='steelblue')
                ylim = ax4_SCZ[r_id_SCZ,c_id_SCZ].get_ylim()
                ax4_SCZ[r_id_SCZ,c_id_SCZ].plot([TH_min_counts, TH_min_counts], ylim,'-k')
            elif opt_style == 'dots':
                ax4_SCZ[r_id_SCZ,c_id_SCZ].scatter(range(1,len(CPC[d_id])+1),np.sort(CPC[d_id]),color = 'steelblue',s=3)
                xlim = ax4_SCZ[r_id_SCZ,c_id_SCZ].get_xlim()
                ax4_SCZ[r_id_SCZ,c_id_SCZ].plot(xlim,[TH_min_counts, TH_min_counts], '-k')
            ax4_SCZ[r_id_SCZ,c_id_SCZ].title.set_text(donor+' ('+disease[d_id]+')')
            
    for g_id in range(0,len(G)): 
        d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
        if opt_style == 'histogram':
            ax4[g_id].hist(np.hstack(CPC[d_IDs]),bins=100,density=True,color='steelblue')
            ylim = ax4[g_id].get_ylim()
            ax4[g_id].plot([TH_min_counts, TH_min_counts], ylim,'-r')
        elif opt_style == 'dots':
            colors_def = cm.rainbow(np.linspace(0, 1, len(d_IDs)))
            for i,d in enumerate(d_IDs):
                ax4[g_id].scatter(range(1,len(np.hstack(CPC[d_IDs[i]]))+1),np.sort(np.hstack(CPC[d_IDs[i]])),color = colors_def[i],s=5, label=sample_IDs[d])
            xlim = ax4[g_id].get_xlim()
            ax4[g_id].plot(xlim,[TH_min_counts, TH_min_counts], '-k')
            ax4[g_id].legend(numpoints=1)
        ax4[g_id].title.set_text('All '+G[g_id]+' donors')
    # if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
    #     for r_id_del in range(r_id_CTRL,n_rows_CTRL):
    #         for c_id_del in range(c_id_CTRL,n_cols):
    #             f4_CTRL.delaxes(ax4_CTRL[r_id_del,c_id_del])
    # if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
    #     for r_id_del in range(r_id_SCZ,n_rows_SCZ):
    #         for c_id_del in range(c_id_SCZ+1,n_cols):
    #             f4_SCZ.delaxes(ax4_SCZ[r_id_del,c_id_del])
    if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
        for r_id_del in range(r_id_CTRL,n_rows_CTRL):
            for c_id_del in range(c_id_CTRL+1,n_cols):
                f4_CTRL.delaxes(ax4_CTRL[r_id_del,c_id_del])
    if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
        for r_id_del in range(r_id_SCZ,n_rows_SCZ):
            for c_id_del in range(c_id_SCZ+1,n_cols):
                f4_SCZ.delaxes(ax4_SCZ[r_id_del,c_id_del])
    if opt_style == 'histogram':        
        f4_CTRL.text(0.5, 0.04, layer_info + ' count depth ('+layer_info+' counts per cell ID)', ha='center')
        f4_SCZ.text(0.5, 0.04, layer_info + ' count depth ('+layer_info+' counts per cell ID)', ha='center')
        f4.text(0.5, 0.04, layer_info + ' count depth ('+layer_info+' counts per cell ID)', ha='center')
        f4_CTRL.text(0.01, 0.5, 'frequency', va='center', rotation='vertical')
        f4_SCZ.text(0.01, 0.5, 'frequency', va='center', rotation='vertical')
        f4.text(0.01, 0.5, 'frequency', va='center', rotation='vertical')
        f4_CTRL.suptitle(layer_info+' count depth distribution',fontsize = FS+5)
        f4_SCZ.suptitle(layer_info+' count depth distribution',fontsize = FS+5)
        f4.suptitle(layer_info+' count depth distribution',fontsize = FS+5)
    elif opt_style == 'dots':
        f4_CTRL.text(0.5, 0.04, 'Cells', ha='center')
        f4_SCZ.text(0.5, 0.04, 'Cells', ha='center')
        f4.text(0.5, 0.04, 'Cells', ha='center')
        f4_CTRL.text(0.01, 0.5,layer_info + ' count Depth ('+layer_info+' counts per cell ID)', va='center', rotation='vertical')
        f4_SCZ.text(0.01, 0.5,layer_info + ' count Depth ('+layer_info+' counts per cell ID)', va='center', rotation='vertical')
        f4.text(0.01, 0.5,layer_info + ' count Depth ('+layer_info+' counts per cell ID)', va='center', rotation='vertical')
        f4_CTRL.suptitle(layer_info + ' cells vs. count depth',fontsize = FS+5)
        f4_SCZ.suptitle(layer_info + ' cells vs. count depth',fontsize = FS+5)
        f4.suptitle(layer_info + ' cells vs. count depth',fontsize = FS+5)
    if opt_save:
        f4_CTRL.savefig(path_results+'count_depth_'+opt_style+'_'+layer_info+'_all_CTRL_samples.pdf', bbox_inches='tight')
        f4_SCZ.savefig(path_results+'count_depth_'+opt_style+'_'+layer_info+'_all_SCZ_samples.pdf', bbox_inches='tight')
        f4.savefig(path_results+'count_depth_'+opt_style+'_'+layer_info+'_all_samples.pdf', bbox_inches='tight')
    plt.close(f4_CTRL)
    plt.close(f4_SCZ)
    plt.close(f4)
    plt.clf()

def plot_number_cells_per_sample(TH_min_cells_per_sample,path_results,filtering_str):

    os.chdir(path_results)
    df = pd.read_pickle("df_"+filtering_str+"_result"+add_str)
    plt = loadPltSettings(20,10)

    fig, axes = plt.subplots(2,1,figsize=(8,8),sharey=True,sharex=True)
    fig.subplots_adjust(hspace = .05, wspace = 0.05)  

    df[df['Disease']=='CTRL'].hist(x='Donor ID', y='Number of Cells',color='steelblue',bins=15,ax=axes[0], edgecolor='none')
    df[df['Disease']=='SCZ'].hist(x='Donor ID', y='Number of Cells',color='darkorange',bins=15,ax=axes[1], edgecolor='none')
    
    for ax in axes: 
        ax.set_ylabel('Number of cells')
        ax.set_title('')
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_xlabel('Donor ID')
        ylim = ax.get_ylim()
        ax.plot([TH_min_cells_per_sample, TH_min_cells_per_sample], ylim,'-b')
    fig.savefig(path_results+'Number_of_cells_after_'+filtering_str+'.pdf', bbox_inches='tight')
    plt.close(fig)
    plt.clf()

def plot_dist_mitoCPC(sample_IDs, disease, PM, RPC, TH_max_fraction_mito_counts, TH_min_counts, path_results, opt_style, opt_save):
    x_max=0
    if disease.count('SCZ')>7:
        n_cols = 7
    else:
        n_cols = np.max([1,disease.count('SCZ')-2])
    n_rows_SCZ = int(np.ceil(np.sum([d=='SCZ' for d in disease])/n_cols))
    n_rows_CTRL = int(np.ceil(np.sum([d=='CTRL' for d in disease])/n_cols))
    plt = loadPltSettings(FS,10)
    f5_CTRL, ax5_CTRL = plt.subplots(n_rows_CTRL, n_cols, figsize=(35,30))
    f5_SCZ, ax5_SCZ = plt.subplots(n_rows_SCZ, n_cols, figsize=(35,30))
    f5, ax5 = plt.subplots(1, 2, figsize=(35,12))
    if opt_style == 'histogram':
        f5_CTRL.subplots_adjust(hspace = .35, wspace = 0.3)  
        f5_SCZ.subplots_adjust(hspace = .35, wspace = 0.3)  
        f5.subplots_adjust(hspace = .35, wspace = 0.1) 
    else:
        f5_CTRL.subplots_adjust(hspace = .55, wspace = 0.15)  
        f5_SCZ.subplots_adjust(hspace = .55, wspace = 0.15)  
        f5.subplots_adjust(hspace = .35, wspace = 0.15)  
    c_id_CTRL=-1
    c_id_SCZ=-1
    r_id_CTRL=0
    r_id_SCZ=0
    for d_id,donor in enumerate(sample_IDs):
        if disease[d_id] =='CTRL':
            c_id_CTRL = c_id_CTRL + 1
            if ((c_id_CTRL)%n_cols==0) & (c_id_CTRL > 0):
                r_id_CTRL = r_id_CTRL + 1
                c_id_CTRL = 0
            if opt_style == 'histogram':
                ax5_CTRL[r_id_CTRL,c_id_CTRL].hist(PM[d_id],bins=100,density=True,color='steelblue')
                ylim = ax5_CTRL[r_id_CTRL,c_id_CTRL].get_ylim()
                if TH_max_fraction_mito_counts<1:
                    ax5_CTRL[r_id_CTRL,c_id_CTRL].plot([TH_max_fraction_mito_counts,TH_max_fraction_mito_counts], ylim,'-r')
                x_max = np.max([x_max,np.max(PM[d_id])])
                ax5_CTRL[r_id_CTRL,c_id_CTRL].set_xlim([0,x_max])
                ax5_CTRL[r_id_CTRL,c_id_CTRL].ticklabel_format(axis = 'x',style='plain')
                #how large is the fraction of cells having percentage of mitochondrial counts > TH_max_fraction_mito_counts*100%??
                if TH_max_fraction_mito_counts<1:
                    ax5_CTRL[r_id_CTRL,c_id_CTRL].text(10,0.9*ylim[1],'Cells with PM>'+str(TH_max_fraction_mito_counts*100)+'%:',fontsize=17)
                    ax5_CTRL[r_id_CTRL,c_id_CTRL].text(50,0.8*ylim[1],str(np.round(100*np.sum(PM[d_id]>TH_max_fraction_mito_counts*100)/len(PM[d_id]),2))+'%',fontsize=17)
            elif opt_style == 'dots':
                RPC_sorted, PM_sorted = map(list,zip(*sorted(zip(RPC[d_id], PM[d_id]))))
                ax5_CTRL[r_id_CTRL,c_id_CTRL].scatter(RPC_sorted,PM_sorted,color='steelblue',s=3)
                xlim = ax5_CTRL[r_id_CTRL,c_id_CTRL].get_xlim()
                ax5_CTRL[r_id_CTRL,c_id_CTRL].plot(xlim,[TH_min_counts, TH_min_counts], '-r', linewidth=1)
            ax5_CTRL[r_id_CTRL,c_id_CTRL].title.set_text(donor+' ('+disease[d_id]+')')
        else:
            c_id_SCZ = c_id_SCZ + 1
            if ((c_id_SCZ)%n_cols==0) & (c_id_SCZ > 0):
                r_id_SCZ = r_id_SCZ + 1
                c_id_SCZ = 0
            if opt_style == 'histogram':
                ax5_SCZ[r_id_SCZ,c_id_SCZ].hist(PM[d_id],bins=100,density=True,color='steelblue')
                ylim = ax5_SCZ[r_id_SCZ,c_id_SCZ].get_ylim()
                if TH_max_fraction_mito_counts<1:
                    ax5_SCZ[r_id_SCZ,c_id_SCZ].plot([TH_max_fraction_mito_counts, TH_max_fraction_mito_counts], ylim,'-r')
                ax5_SCZ[r_id_SCZ,c_id_SCZ].plot([5,5], ylim,'-r')
                x_max = np.max([x_max,np.max(PM[d_id])])
                ax5_SCZ[r_id_SCZ,c_id_SCZ].set_xlim([0,x_max])
                ax5_SCZ[r_id_SCZ,c_id_SCZ].ticklabel_format(axis = 'x',style='plain')
                #how large is the fraction of cells having percentage of mitochondrial counts > TH_max_fraction_mito_counts*100%??
                if TH_max_fraction_mito_counts<1:
                    ax5_SCZ[r_id_SCZ,c_id_SCZ].text(10,0.9*ylim[1],'Cells with PM>'+str(TH_max_fraction_mito_counts*100)+'%:',fontsize=17)
                    ax5_SCZ[r_id_SCZ,c_id_SCZ].text(50,0.8*ylim[1],str(np.round(100*np.sum(PM[d_id]>TH_max_fraction_mito_counts*100)/len(PM[d_id]),2))+'%',fontsize=17)
            elif opt_style == 'dots':
                RPC_sorted, PM_sorted = map(list,zip(*sorted(zip(RPC[d_id], PM[d_id]))))
                ax5_SCZ[r_id_SCZ,c_id_SCZ].scatter(RPC_sorted,PM_sorted,color='steelblue',s=3)
                xlim = ax5_SCZ[r_id_SCZ,c_id_SCZ].get_xlim()
                ax5_SCZ[r_id_SCZ,c_id_SCZ].plot(xlim,[TH_min_counts, TH_min_counts], '-r', linewidth=1)
                
            ax5_SCZ[r_id_SCZ,c_id_SCZ].title.set_text(donor+' ('+disease[d_id]+')')
    for g_id in range(0,len(G)): 
        d_IDs = [i for i in range(0,len(disease)) if disease[i]==G[g_id]]
        if opt_style == 'histogram':
            ax5[g_id].hist(np.hstack(PM[d_IDs]),bins=100,density=True,color='steelblue')
            ylim = ax5[g_id].get_ylim()
            if TH_max_fraction_mito_counts<1:
                ax5[g_id].plot([TH_max_fraction_mito_counts,TH_max_fraction_mito_counts], ylim,'-r')
            x_max = np.max([x_max,np.max(np.hstack(PM[d_IDs]))])
            ax5[g_id].set_xlim([0,x_max])
            ax5[g_id].ticklabel_format(axis = 'x',style='plain')
            if TH_max_fraction_mito_counts<1:
                #how large is the fraction of cells having percentage of mitochondrial counts > TH_max_fraction_mito_counts*100%??
                ax5[g_id].text(10,0.9*ylim[1],'Cells with PM>'+str(TH_max_fraction_mito_counts*100)+'%:',fontsize=17)
                ax5[g_id].text(50,0.8*ylim[1],str(np.round(100*np.sum(np.hstack(PM[d_IDs])>TH_max_fraction_mito_counts*100)/len(np.hstack(PM[d_IDs])),2))+'%',fontsize=17)
        elif opt_style == 'dots':
            RPC_sorted, PM_sorted = map(list,zip(*sorted(zip(np.hstack(RPC[d_IDs]), np.hstack(PM[d_IDs])))))
            ax5[g_id].scatter(RPC_sorted,PM_sorted,color='steelblue',s=3)
            xlim = ax5[g_id].get_xlim()
            ax5[g_id].plot(xlim,[TH_min_counts, TH_min_counts], '-r', linewidth=1)
        ax5[g_id].title.set_text('All ' + G[g_id] + ' donors')
    
    # if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
    #     for r_id_del in range(r_id_CTRL,n_rows_CTRL):
    #         for c_id_del in range(c_id_CTRL,n_cols):
    #             f5_CTRL.delaxes(ax5_CTRL[r_id_del,c_id_del])
    # if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
    #     for r_id_del in range(r_id_SCZ,n_rows_SCZ):
    #         for c_id_del in range(c_id_SCZ+1,n_cols):
    #             f5_SCZ.delaxes(ax5_SCZ[r_id_del,c_id_del])
    if ((r_id_CTRL)*n_cols+c_id_CTRL+1) < (n_rows_CTRL*n_cols):
        for r_id_del in range(r_id_CTRL,n_rows_CTRL):
            for c_id_del in range(c_id_CTRL+1,n_cols):
                f5_CTRL.delaxes(ax5_CTRL[r_id_del,c_id_del])
    if ((r_id_SCZ)*n_cols+c_id_SCZ+1) < (n_cols*n_rows_SCZ):
        for r_id_del in range(r_id_SCZ,n_rows_SCZ):
            for c_id_del in range(c_id_SCZ+1,n_cols):
                f5_SCZ.delaxes(ax5_SCZ[r_id_del,c_id_del])
    f5_CTRL.text(0.5, 0.04, 'Percentage of mitochondrial transcripts per cell (PM)', ha='center')
    f5_SCZ.text(0.5, 0.04, 'Percentage of mitochondrial transcripts per cell (PM)', ha='center')
    f5.text(0.5, 0.01, 'Percentage of mitochondrial transcripts per cell (PM)', ha='center')
    if opt_style == 'histogram':
        f5_CTRL.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
        f5_SCZ.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
        f5.text(0.01, 0.5, 'Frequency', va='center', rotation='vertical')
        f5_CTRL.suptitle('Mitochondrial transcript distribution',fontsize = FS+5)
        f5_SCZ.suptitle('Mitochondrial transcript distribution',fontsize = FS+5)
        f5.suptitle('Mitochondrial transcript distribution',fontsize = FS+5)
    elif opt_style == 'dots':
        f5_CTRL.text(0.01, 0.5,'Count depth', va='center', rotation='vertical')
        f5_SCZ.text(0.01, 0.5,'Count Depth', va='center', rotation='vertical')
        f5.text(0.01, 0.5,'Count Depth', va='center', rotation='vertical')
        f5_CTRL.suptitle('Mitochondrial content',fontsize = FS+5)
        f5_SCZ.suptitle('Mitochondrial content',fontsize = FS+5)
        f5.suptitle('Mitochondrial content',fontsize = FS+5)
    if opt_save:
        f5_CTRL.savefig(path_results+'distribution_percentage_mito_counts_'+opt_style+'_all_CTRL_samples.pdf', bbox_inches='tight')
        f5_SCZ.savefig(path_results+'distribution_percentage_mito_counts_'+opt_style+'_all_SCZ_samples.pdf', bbox_inches='tight')
        f5.savefig(path_results+'distribution_percentage_mito_counts_'+opt_style+'_all_samples.pdf', bbox_inches='tight')
    plt.close(f5)
    plt.close(f5_SCZ)
    plt.close(f5_CTRL)
    plt.clf()

def plot_data(PD, PD_updated, TH_min_number_counts, TH_min_number_s_counts, TH_min_number_us_counts, TH_max_fraction_mito_counts, TH_min_number_genes, TH_min_cells_per_sample, path_project,  path_summary_metrics, opt_TH_filtered, opt_doublets_removed, opt_manual_doublet_TH_adjustment, opt_stats_data, opt_save, opt_D_version, opt_pagoda_filtering):

    filtering_str = get_filtering_str(opt_TH_filtered,opt_doublets_removed,opt_manual_doublet_TH_adjustment,opt_D_version,opt_pagoda_filtering)

    path_data_raw, path_code, path_results_raw = get_paths(path_project,opt_pagoda_filtering, 'raw',0)
    path_data_filt, path_code, path_results_filt = get_paths(path_project,opt_pagoda_filtering, 'preprocessed',0)
    if filtering_str == '':
        D = get_raw_data([],[],False,path_project,path_code,path_data_raw,path_summary_metrics)
        path_results = path_results_raw
    elif filtering_str == '_TH_filtered':
        D, _, _, _ = get_TH_filtered_data([], [], path_code, path_data_raw, path_data_filt, path_results, False,True,opt_pagoda_filtering)
        path_results = path_results_filt
    else:
        D, _, _, _, _, _ = get_TH_and_D_filtered_data([], [], path_code, path_data_raw, path_data_filt, path_results, opt_pagoda_filtering, False, True, True, opt_manual_doublet_TH_adjustment, False, TH_min_cells_per_sample, opt_D_version)
        path_results = path_results_filt
    MCG, CPG, min_CPG, GPC, CPC, sCPC, usCPC, MCPC, PMito, RCPC, PRibo, mean_counts_per_cell, std_counts_per_cell, RLD, med_gpc, df_stats = get_stats_per_sample(D,path_results,path_code,opt_stats_data,filtering_str,opt_save)
    
    sample_IDs,_ ,disease, sex = get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
    
    D.close() 

    layer_info = [""]

    ## plot count depth vs. genes per cell ID with no or gradual/ thrsholded mitochondrial counts per cell color code
    for l_id,layer in enumerate(layer_info):
        if layer=="spliced":
            cpc = sCPC
        elif layer=="unspliced":
            cpc = usCPC
        else:
            cpc = CPC
            opt_colorCode = 'sample'
            plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, [], TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode, opt_save)

        opt_colorCode = 'group' 
        plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, [], TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode, opt_save)
        
        opt_colorCode = 'none' 
        plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, [], TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode, opt_save)

        opt_colorCode = 'mitoRPC_gradient'
        plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, [], TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode, opt_save) 
        
        if TH_max_fraction_mito_counts<1:
            opt_colorCode = 'mitoRPC_cutoff'
            plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, [], TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode,  opt_save)
        
        #To do: needs to be fixed within plot_CPC_vs_GPC
        if opt_doublets_removed==True:
            opt_colorCode = 'doublet_prediction'
            if opt_manual_doublet_TH_adjustment:
                plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, PD_updated, TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode,  opt_save)
            else:
                plot_CPC_vs_GPC(sample_IDs, disease, cpc, GPC, PMito, PD, TH_min_number_counts, TH_min_number_genes, TH_max_fraction_mito_counts, path_results, layer, opt_colorCode, opt_save)
    
    ### plot histogram of mean counts per gene
    plot_MCG(sample_IDs, disease, MCG, path_results, opt_save)
                            
    ## plot distribution of genes per cell ID:
    ## - low number of genes per barcode: quiescent cells
    ## - high number of genes per barcode: doublets
    plot_GPC(sample_IDs, disease, sex, GPC, TH_min_number_genes, path_results, opt_save)
    
    ## plot count depth: distribution of counts per CellID (barcode) for each sample
    ## - low count per barcode: quiescenct cells
    ## - high count per barcode: doublets and large cells (also for nuclei?)
    for l_id,layer in enumerate(layer_info):
        if layer=="spliced":
            cpc = sCPC
        elif layer=="unspliced":
            cpc = usCPC
        else:
            cpc = CPC
        opt_style = 'histogram'
        plot_CPC(sample_IDs, disease, cpc, path_results, TH_min_number_counts,layer, opt_style, opt_save)
        opt_style = 'dots'
        plot_CPC(sample_IDs, disease, cpc, path_results, TH_min_number_counts, layer, opt_style, opt_save)

    ## plot distribution of the percentage of mitochondrial counts per cell
    ## - cells with high percentage of mitochondrial counts can be in respiratory process 
    ##   or belong to cells whose membranes are broken such that mRNA has leaked out and 
    ##   only mitochondrial mRNA is still conserved
    opt_style = 'histogram'
    plot_dist_mitoCPC(sample_IDs, disease, PMito, CPC, TH_max_fraction_mito_counts, TH_min_number_counts, path_results, opt_style, opt_save)
    opt_style = 'dots'
    plot_dist_mitoCPC(sample_IDs, disease, PMito, CPC, TH_max_fraction_mito_counts, TH_min_number_counts, path_results, opt_style, opt_save)   
    
    #plot histogram of mito counts per cell for each sample
    #To do:
    #plot_mitoCPC(sample_IDs, disease,MCPC,PMito,TH_max_fraction_mito_counts,path_results,opt_save)
    
    #plot histogram of number of cells per sample after filtering:
    if path_results != path_results_raw:
        plot_number_cells_per_sample(TH_min_cells_per_sample,path_results,filtering_str)
    
    gc.collect()
    #optional: count depth vs. protein coding fraction, count depth vs. ribosomal fraction 

def plot_doublets_scrublet(scrub,donor,path_results):  
    plt = loadPltSettings(12,10)
    scrub.plot_histogram()
    plt.savefig(path_results+'/doublets_scrublet_hist_'+donor+'.pdf')
    plt.close()
    plt.clf()
    #manual adjustment of threshold: if automatic threshold detection doesn't work well, you can adjust the threshold with the call_doublets() function. For example:
    #scrub.call_doublets(threshold=0.25)
    #plot doublet predictions on 2D embedding --> predicted doublets should co-localize in dinstinct states
    plt = loadPltSettings(12,10)
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP',order_points=True)
    plt.savefig(path_results+'/doublets_scrublet_UMAP_'+donor+'.pdf')
    plt.close()
    plt.clf()

def plot_distributions(df_sel,path_results,var):
    plt = loadPltSettings(4,3)
    if var in ["age","PMI_h"]:
        sns.kdeplot(data=df_sel.loc[df_sel.scz_status=="ctrl"],x=var,color="steelblue",label = "CTRL", linewidth = 3)
        lm = sns.kdeplot(data=df_sel.loc[df_sel.scz_status=="scz"],x=var,color="darkorange",label = "SCZ", linewidth = 3)
        lm.axes = remove_frame_keep_axes(lm.axes)
        lm.axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        def get_abundance(df_sel,s,d):
            if d=="CTRL":
                num = np.sum(df_sel[df_sel.scz_status=="ctrl"]["gender"]==s)
            elif d=="SCZ":
                num = np.sum(df_sel[df_sel.scz_status=="scz"]["gender"]==s)
            return num
        sex = ["m","m","f","f"]
        disease = ["CTRL","SCZ","CTRL","SCZ"]
        abu = [get_abundance(df_sel,s,d) for s,d in zip(sex,disease)]
        df = pd.DataFrame({"disease":disease,"abundance":abu,"sex":sex})
        ax = plt.bar(data=df,height="abundance", x="sex",color=["steelblue","darkorange","steelblue","darkorange"])
    if not os.path.exists(path_results):
        os.makedirs(path_results)         
    plt.savefig(path_results+var+'_distribution_per_group.pdf', bbox_inches='tight')
    plt.close()
    plt.clf() 

def plot_percentage_cells_discarded(file_str,path_results,opt_together):
    os.chdir(path_results )
    df = pd.read_pickle("df_"+file_str+"_result"+add_str)
    plt = loadPltSettings(20,10)
    if opt_together:
        fig, ax = plt.subplots()
        a_heights, a_bins = np.histogram(df[df['Disease']=='CTRL']['Percentage Cells discarded'])
        b_heights, b_bins = np.histogram(df[df['Disease']=='SCZ']['Percentage Cells discarded'], bins=a_bins)
        width = (a_bins[1] - a_bins[0])/3
        ax.bar(a_bins[:-1], a_heights, width=width, facecolor='steelblue')
        ax.bar(b_bins[:-1]+width, b_heights, width=width, facecolor='darkorange')
        ax.set_ylabel('Absolute frequency')
        ax.set_title('')
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_xlabel('Percentage cells discarded')
    else:
        fig, axes = plt.subplots(2,1,figsize=(8,8),sharey=True,sharex=True)
        fig.subplots_adjust(hspace = .05, wspace = 0.05)  
        #df['Percentage Cells discarded'].hist(bins=20,by=df['Disease'],color={'blue','red'},ax=axes)
        df[df['Disease']=='CTRL'].hist('Percentage Cells discarded',color='steelblue',bins=15,ax=axes[0], edgecolor='none')
        df[df['Disease']=='SCZ'].hist('Percentage Cells discarded',color='darkorange',bins=15,ax=axes[1], edgecolor='none')
        for ax in axes: 
            ax.set_ylabel('Absolute frequency')
            ax.set_title('')
            ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_xlabel('Percentage cells discarded')
    #fig.savefig(path_results+'Percentage_cells_discrded_after_'+file_str+'.png', bbox_inches='tight')
    fig.savefig(path_results+'Percentage_cells_discarded_after_'+file_str+'.pdf', bbox_inches='tight')
    plt.close(fig)
    plt.clf()

def plot_UMAP_removed_cells(path_figure,path_loom_file,loom_file_name,opt_pagoda_filtering):
    # Create an array with the colors you want to use
    if opt_pagoda_filtering:
        my_filtering_colors = ["#003366","#CC0066","#078383", "#939393"] #P filtered', 'TH filtered', 'D filtered','not filtered'
    else:
        my_filtering_colors = ["#CC0066","#078383","#939393"] #'TH filtered', 'D filtered','not filtered'

    os.chdir(path_loom_file)

    aD = anndata.read_loom(loom_file_name, X_name='matrix') 
    aD.var_names_make_unique()
    
    #remove x and y chromosome information:
    DF_xy_chr = pd.read_csv(path_loom_file+'XYchr_genes.csv',delimiter=',',index_col=None)
    non_xy_chr_genes_list =[g for a,g in zip(aD.var['Accession'],aD.var.index) if a not in DF_xy_chr['Accession'].tolist()]
    aD = aD[:,non_xy_chr_genes_list]  
    
    #remove samples removed:
    bool_cells_kept = np.logical_and(aD.obs['Donor']!='S78', aD.obs['Donor']!='S66')
    aD = aD[bool_cells_kept,:]
    
    #get colors for anndata object metavariables
    aD = get_color_maps_for_anndata(aD)
    
    ##basic normalization
    #normalize each cell by total counts over all genes, so that every cell has the same total count after normalization 
    #normalize to 10 000 reads per cell --> counts become comparable amongst cells
    sc.pp.normalize_total(aD, layers = 'all')
    ##log-transformation of data
    #Computes X=log(X+1), where log denotes the natural logarithm unless a different base is given.
    #chunked=True saves memory
    sc.pp.log1p(aD,chunked=True)#
    #determine HVGs --> changes neighborhood graph calculation
    sc.pp.highly_variable_genes(aD,flavor='cell_ranger',n_top_genes=2000,inplace=True)
    ##embedding the neihborhood graph in 2 dimensions
    ##plot umap (all samples, all samples healthy, all samples SCZ on top of each other) with gene expression of highly variable gene as color code for top 20 highly variable genes
    os.chdir(path_figure)
    #pca
    sc.tl.pca(aD, svd_solver='arpack')
    sc.pl.pca(aD)
    variance_explained_PCA = aD.uns["pca"]["variance_ratio"]
    np.save('PCA_var_explained_'+loom_file_name[0:-5],variance_explained_PCA,allow_pickle=True)
    sc.pp.neighbors(aD, n_neighbors=10, n_pcs=30)
    sc.tl.umap(aD,random_state=192785)
    if opt_pagoda_filtering:
        aD.obs['filtering_status'].cat.reorder_categories(['P filtered', 'TH filtered', 'D filtered','not filtered'], inplace=True)
    else:
        aD.obs['filtering_status'].cat.reorder_categories(['TH filtered', 'D filtered','not filtered'], inplace=True)
    sc.pl.umap(aD, color=['filtering_status'],save="Samples_UMAP_filtering_status.pdf",palette=my_filtering_colors)

def plot_expression_PCA(loom_file_str,file_str1,file_str2,opt_loom,opt_violin_plots,opt_spearman_corr_with_QC_metrics,path_data,path_results_filt):
    os.chdir(path_data+'filtered_loom_formatted_data/')  
    aD_i = anndata.read_loom('Samples'+add_str+loom_file_str, X_name='matrix') 
    figure_str = add_str+'_TH_and_D_adj_filtered'
    aD_i.var_names_make_unique()
    if opt_loom:
        non_xy_chr_genes_list = aD_i.var_names[np.logical_and(aD_i.var['Chromosome']!='Y',aD_i.var['Chromosome']!='X')]
        aD = aD_i[:,non_xy_chr_genes_list]
    ##basic normalization
    #normalize each cell by total counts over all genes, so that every cell has the same total count after normalization 
    sc.pp.normalize_total(aD, layers = 'all')
    #Computes X=log(X+1), where log denotes the natural logarithm unless a different base is given.
    #chunked=True saves memory
    sc.pp.log1p(aD,chunked=True)#
    
    #each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization.
    # PCA of expression values
    sc.tl.pca(aD, n_comps = 30, svd_solver='arpack')
    # color code by diesease, individual & library
    #sc.pl.pca(aD,color=['Disease','Donor'],save='_subset'+add_str+'_TH_and_D_adj_filtered')
    sc.pl.pca(aD,color=['Donor'],save='_donor'+figure_str)
    sc.pl.pca(aD,color=['Library','Sex'],save='_library_sex'+figure_str)
    sc.pl.pca(aD,color=['Disease','PMI_h'],save='_diesese_PMI'+figure_str)
    sc.pl.pca_variance_ratio(aD, log=False, save=figure_str)
    sc.pl.pca_loadings(aD, components=None, show=True, save=figure_str)
    #aD.obsm['X_pca'] #PCA representation of data. (46362, 50)
    #aD.varm['PCs']# The principal components containing the loadings. (30095, 50)
    #aD.uns['pca']['variance_ratio'] # Ratio of explained variance. 
    #aD.uns['pca']['variance']#Explained variance, equivalent to the eigenvalues of the covariance matrix.
    #To Do: plot PCs seperately e.g. violin plots
    
    #calculate spearman correlation of expression values with quality metrics
    for i in range(0,np.shape(aD.obsm['X_pca'])[1]):
        if i==0:
            d = {'donor_ID_python': aD.obs['Donor'], 'PC1': aD.obsm['X_pca'][:,0]}
            df_ePCA = pd.DataFrame(data=d)
        else:
            df_ePCA['PC'+str(i+1)] = aD.obsm['X_pca'][:,i]
    os.chdir(path_data)
    
    df_qc_metrics = pd.read_csv('metrics_summary_tidy.csv')
    #integrate % removed cells during TH filtering per sample
    #integrate % removed doublets per sample
    os.chdir(path_results_filt)
    df_stats_TH = pd.read_pickle('df_'+file_str1+'_result'+add_str)
    df_stats_D = pd.read_pickle('df_'+file_str2+'_result'+add_str)
    df_stats_TH.rename(columns = {'Donor ID':'donor_ID_python', 'Percentage Cells discarded':'P_cells_del_TH_filt'},inplace=True)
    df_stats_D.rename(columns = {'Donor ID':'donor_ID_python', 'Percentage Cells discarded':'P_cells_del_D_filt'},inplace=True)
    df_qc_metrics = pd.merge(df_qc_metrics,df_stats_TH[['P_cells_del_TH_filt','donor_ID_python']],on='donor_ID_python')
    df_qc_metrics = pd.merge(df_qc_metrics,df_stats_D[['P_cells_del_D_filt','donor_ID_python']],on='donor_ID_python')
    df_qc_and_pca = pd.merge(df_ePCA,df_qc_metrics,on='donor_ID_python')
    df_stats_TH.rename(columns = {'donor_ID_python':'Donor'},inplace=True)
    df_stats_D.rename(columns = {'donor_ID_python':'Donor'},inplace=True)
    aD.obs = pd.merge(aD.obs,df_stats_TH[['P_cells_del_TH_filt','Donor']],on='Donor')
    aD.obs = pd.merge(aD.obs,df_stats_D[['P_cells_del_D_filt','Donor']],on='Donor')
    
    if opt_violin_plots:
        #plot PCs seperately as bee swarm / violin plots
        PCs_str = ["PC1","PC2","PC3","PC4"]
        color_groups = [["gender","Library"],["P_cells_del_TH_filt","P_cells_del_D_filt"]]
        sns.set_context("paper",font_scale=3)
        #sns.set_context("paper", rc={"font.size":30,"axes.titlesize":30,"axes.labelsize":30}) 
        for cg_id1 in range(np.shape(color_groups)[0]):
            f8, ax8 = plt.subplots(len(PCs_str), np.shape(color_groups)[1], sharex=True, sharey=True, figsize=(35,35))
            for pc_id,pc in enumerate(PCs_str):
                for cg_id2 in range(np.shape(color_groups)[1]):
                    g=sns.violinplot(x="Group", y=pc,  hue=color_groups[cg_id1][cg_id2], data=df_qc_and_pca, dodge = True, ax=ax8[pc_id,cg_id2], palette="Set2")
                    #g=sns.swarmplot(x="Group", y=pc,  hue=cg, data=df_qc_and_pca, dodge = True, ax=ax8[pc_id,cg_id], size=1, palette="Set2", edgecolor = 'gray')
                    g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            f8.savefig(path_results_filt+'/figures/pcs_violin_beeswarm'+figure_str+'.png', bbox_inches='tight')
            plt.close(f8)
            plt.clf()
    
    if opt_spearman_corr_with_QC_metrics:
        Corr_spearman, qc_ids = calculate_spearman_corr_with_QC_metrics(df_qc_and_pca)
        
        # grouped barplot:
        barWidth = 0.06
        n_pcs_to_plot = 5
        colors = cm.tab20(range(np.max(np.shape(qc_ids))))
        for qc_id in range(0,np.max(np.shape(qc_ids))):
            # Set position of bar on X axis
            if qc_id==0:
                r = np.arange(len(Corr_spearman[0,0:n_pcs_to_plot]))
            else:
                r = [x + barWidth for x in r]
            plt.bar(r, Corr_spearman[qc_id,0:n_pcs_to_plot], color=colors[qc_id], width=barWidth, edgecolor='white', label=df_qc_and_pca.columns[qc_ids][qc_id])
        # Add xticks on the middle of the group bars
        plt.xlabel('Expression value principle components', fontweight='bold')
        plt.ylabel('Spearman correlation', fontweight='bold')
        xticks_str = list()
        for pc_id in range(0,n_pcs_to_plot):
            xticks_str.append('PC'+str(pc_id+1)+' \n('+str(np.round((aD.uns['pca']['variance_ratio'][pc_id]*100),2))+'%)')
        plt.xticks([r + barWidth for r in range(len(Corr_spearman[0,0:n_pcs_to_plot]))],xticks_str)
        # Create legend & Show graphic
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        #plt.show() 
        plt.savefig(path_results_filt+'figures/'+'SpearmanCorr_qc_metrics_with_exprPCs_'+figure_str+'.png', bbox_inches='tight')
        plt.savefig(path_results_filt+'figures/'+'SpearmanCorr_qc_metrics_with_exprPCs_'+figure_str+'.pdf', bbox_inches='tight')
        plt.close()
        plt.clf()


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