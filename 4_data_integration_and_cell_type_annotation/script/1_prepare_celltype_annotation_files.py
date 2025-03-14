# prepare files for cell type annotation based on filtered scRNAseq/ snRNAseq data
# Created on Wed Nov 11 12:48:51 2020
# @author: Lisa Bast
# version: 0.0.3

import os, sys
import loompy
import random
import pandas as pd 
import numpy as np

random.seed(30)

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("4_data_integration_and_cell_type_annotation\\script","")

path_data, _, path_results = ut.get_paths(path_project,"Cell type annotation")

# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(len(sys.path)+1, path_code)

## specify settings:
ref_data_set='Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq' 
ref_data_set_short = 'ABM_MCA'

opt_pagoda_filtering=False
opt_save_sample_as_csv = False 
if opt_save_sample_as_csv:
    opt_D_version = 'v1'
opt_format_refData = True 
opt_whole_query_dataset = True #False 
opt_ref_data_as_query_data=False # if true ref data is used for testing 

if opt_save_sample_as_csv:
    num_cells_in_query=2500
    if opt_whole_query_dataset==False and opt_ref_data_as_query_data == False:
        sample_ID = 'S2'#'S48'
        
if opt_format_refData:
    opt_red_ref_dataset = True #False # 
    if opt_red_ref_dataset:
        #red_method = 'HVG_with_seurat_v3' 
        red_method = 'HVG_with_cell_ranger'
        #red_method = 'highest_dropout_rate'
        n_genes = [5000]#[500,1000,2000,3000]#,5000]
    else:
        red_method = ''
        n_genes = [0]
    opt_num_CT_clusters = [51,76] #18 or 51 or 76 or 120
    opt_all_genes = False #if false top 500 highly variable genes in ref dataset are kept in ref data; 
                          #if true all genes also measured in test data are taken

add_str = "_cellranger"

## save query data set:
samples_selected = []
IDs_cell_sel = []
if opt_save_sample_as_csv:
    #get sanple(s) and store matrix as .csv file
    if opt_ref_data_as_query_data:
        path_query_data = path_data+"reference_data_sets/"+ref_data_set_short
        fileName = "matrix"
    else:
        path_query_data = path_data
        if opt_whole_query_dataset:
            fileName = "Samples"+add_str+"_TH_and_D_adj_"+opt_D_version+"_filtered"
    os.chdir(path_query_data)  
    D = loompy.connect(fileName+".loom")
    if opt_whole_query_dataset:
        if opt_ref_data_as_query_data:
            pd.DataFrame(D[:,:],index=D.ra.Gene).to_csv(path_query_data+"/"+fileName+".csv")#, header=False)
            D.close()
        else:
            sample_IDs,_,_,_ = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
            D.close()
            # make sure S78 and S66 are not included
            sample_IDs = sample_IDs[((sample_IDs!='S78') & (sample_IDs!='S66'))]
            for s_id,sample in enumerate(sample_IDs):
                fileName_i = sample+add_str+"_TH_and_D_adj_"+opt_D_version+"_filtered"
                if opt_D_version=='v1':
                    D = loompy.connect(sample+add_str+"_TH_and_D_adj_filtered.loom")
                else:
                    D = loompy.connect(fileName_i+".loom")
                #store num_cells_in_query cells at a time in a .csv file
                n_files = int(np.ceil(np.shape(D)[1]/num_cells_in_query))
                for f_id in range(1,n_files+1):
                    cID_start = 0+(f_id-1)*num_cells_in_query
                    cID_end = f_id*num_cells_in_query
                    if (np.shape(D)[1]<cID_end):
                        cID_end = np.shape(D)[1]
                    pd.DataFrame(D[:,cID_start:cID_end],index=D.ra.Gene).to_csv(path_query_data+"/"+ fileName_i + "_"+ str(f_id) +".csv")#, header=False)
                D.close()
    else:
        #randomly pick num_cells_in_query cells
        IDs_cell_sel = np.sort(random.sample(range(np.shape(D[:,:])[1]),num_cells_in_query)).flatten()
        #counts_per_cell = np.sum(D[:,:],axis=0)
        #idx = np.sort((-counts_per_cell).argsort()[:num_cells_in_query])
        if opt_ref_data_as_query_data:
            group_str = np.repeat('train',len(D.ca.obs_names))
            group_str[IDs_cell_sel] = 'test'
            D.set_attr("group",group_str,axis=1)
            obs_names = D.ca['obs_names']
            group_names = D.ca['group']
            col_attrs = { "cell_name": obs_names, "group": group_names}
            genenames = D.ra['var_names']
            row_attrs = { "Gene": genenames}
            loompy.create('matrix_split_into_test_and_train.loom', D[:,:], row_attrs, col_attrs)
            df_matrix_test = pd.DataFrame(D[:,IDs_cell_sel],index=D.ra.var_names)
            df_matrix_test.to_csv(path_query_data+"/"+fileName+"_"+str(num_cells_in_query)+"_cells_test.csv")#, header=False)
            IDs_cell_train = np.setdiff1d(range(0,len(D.ca.obs_names)),IDs_cell_sel)
            df_matrix_train = pd.DataFrame(D[:,IDs_cell_train],index=D.ra.var_names)
            df_matrix_train.to_csv(path_query_data+"/"+fileName+"_"+str(len(D.ca.obs_names))+"_cells_train.csv")#, header=False)
            samples_selected = D.ca['obs_names'][IDs_cell_sel]
            D.close()
            #store GT labels:
            os.chdir(path_project+'/data/reference_data_sets/'+ref_data_set_short)
            metadata = pd.read_csv('metadata.csv')
            
            metadata = ut.merge_CT_cluster_of_metadata(ref_data_set,metadata,opt_num_CT_clusters)
            df_map = metadata[['sample', 'cell type']].copy()
            df_map_test = df_map.copy()  
            df_map_test['index_whole_data_set'] = df_map_test.index
            #save df with test set and train set
            df_map_test_sel = df_map_test[df_map_test['sample'].isin(samples_selected)]
            df_map_test_sel.reset_index(inplace=True)
            df_map_test_sel.drop(columns=['index'],inplace=True)
            os.chdir(path_query_data)   
            df_map_test_sel.to_csv(path_query_data+"/"+fileName+"_"+str(num_cells_in_query)+"_cells_GT_test.csv",index=True)
            df_map_train_sel = df_map_test[~df_map_test['sample'].isin(samples_selected)]
            df_map_train_sel.reset_index(inplace=True)
            df_map_train_sel.drop(columns=['index'],inplace=True)
            os.chdir(path_query_data)   
            df_map_train_sel.to_csv(path_query_data+"/"+fileName+"_"+str(np.shape(df_map_train_sel)[0])+"_cells_GT_train.csv",index=True)
        else:
            pd.DataFrame(D[:,IDs_cell_sel],index=D.ra.Gene).to_csv(path_query_data+"/"+fileName+"_"+str(num_cells_in_query)+"_cells.csv")#, header=False)
            D.close()
            
## format and save the reference data set:
if opt_format_refData:
    test_data_path = path_data+"filtered_loom_formatted_data"
    test_data_filename = "Samples"+add_str+"_TH_and_D_adj_filtered.loom"
    for n_g in n_genes:
        ut.format_ref_data(path_project, ref_data_set, ref_data_set_short, opt_num_CT_clusters, opt_all_genes, red_method, test_data_path, test_data_filename, IDs_cell_sel, samples_selected, opt_ref_data_as_query_data, n_g)

        