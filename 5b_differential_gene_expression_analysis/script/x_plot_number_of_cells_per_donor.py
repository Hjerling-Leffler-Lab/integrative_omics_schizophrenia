# about: plot number of cells per donor and cell type
# date: Created on Wed Aug  4 13:42:20 2021
# @author: Lisa Bast
# version: 0.0.2

#to do: fix paths and filenames and delete unccessary stuff

import os
import numpy as np
import loompy
import pandas as pd
import sys
import gc

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("5b_differential_gene_expression_analysis\\script","")

opt_pagoda_filtering = True
opt_plot=False#True#
#opt_celltype_groups = 'scmap'
opt_celltype_groups = 'conos_cluster_based'

opt_df_cells_per_donor_and_CT = 'load' #'create'#
opt_downsample_cells = False#
opt_use_downsampled_loom_files = False
if opt_downsample_cells or opt_use_downsampled_loom_files:
    n_datasets = 10


add_str = '_cellranger'
n_clusters_range = [15,37]

if opt_pagoda_filtering:
    add_str_pagoda = '_pagoda'
else:
    add_str_pagoda = ''

    
for n_clusters in n_clusters_range:
    # read data as anndata object
    path_data, _, path_results = ut.get_paths(path_project,'DEG_ana_preparation') 
    
    path_data_source=path_data
    cl_str = '_'+str(n_clusters)+'_CTs'
    os.chdir(path_data_source)

    str_not_contains = 'downsampled'
    str_start = 'Samples_cellranger'+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups
    str_end = '.loom'
    loom_files = [ f for f in os.listdir( os.curdir ) if os.path.isfile(f) and f.startswith(str_start) and f.endswith(str_end) and str_not_contains not in f]

    if opt_df_cells_per_donor_and_CT == 'create':
        bool_first=True
        for loom_file_name in loom_files:
            #explore files
            #extract number of cells per donor, group and cell type
            D = loompy.connect(loom_file_name)
            CT_str = loom_file_name[len(str_start)+1:-5]
            if bool_first:
                df_donor_disease = pd.DataFrame(data={'Donor':D.ca['Donor'].tolist(),'Disease':D.ca['Disease'].tolist()}).drop_duplicates(ignore_index=True)
                df = pd.DataFrame(data = {'Cell type cluster': CT_str, 'cells_total':np.shape(D)[1], 'cells_SCZ': sum(D.ca['Disease']=='SCZ'), 'cells_CTRL': sum(D.ca['Disease']=='CTRL')},index=[0])
                df_p = pd.DataFrame(data = {'Cell type cluster': CT_str},index=[0])
                for donor in np.unique(D.ca['Donor']).tolist():
                    df[donor]=sum(D.ca['Donor']==donor)
                    df_p[donor]=(sum(D.ca['Donor']==donor)/np.shape(D)[1])*100
                bool_first=False
            else:
                df_add=pd.DataFrame(data = {'Cell type cluster': CT_str, 'cells_total':np.shape(D)[1], 'cells_SCZ': sum(D.ca['Disease']=='SCZ'), 'cells_CTRL': sum(D.ca['Disease']=='CTRL')},index=[0])
                df_p_add = pd.DataFrame(data = {'Cell type cluster': CT_str},index=[0])
                for donor in np.unique(D.ca['Donor']).tolist():
                    df_add[donor]=sum(D.ca['Donor']==donor)
                    df_p_add[donor]=(sum(D.ca['Donor']==donor)/np.shape(D)[1])*100
                df=df.append(df_add,ignore_index=True,sort=False)
                df_p=df_p.append(df_p_add,ignore_index=True,sort=False)
            D.close()
        df = df.set_index(keys='Cell type cluster', drop=True)
        df_p = df_p.set_index(keys='Cell type cluster', drop=True)
        #create dataframe with n_cells per donor
        donors_col = [col for col in df if col.startswith('S')]
        df_donors = df[donors_col].copy().transpose()
        df_donors['Donor'] = df_donors.index
        df_donors=pd.merge(left=df_donors,right=df_donor_disease,on='Donor')
        
        donors_col_p = [col for col in df_p if col.startswith('S')]
        df_p_donors = df_p[donors_col_p].copy().transpose()
        df_p_donors['Donor'] = df_p_donors.index
        df_p_donors=pd.merge(left=df_p_donors,right=df_donor_disease,on='Donor')
        
        #save data frames
        df_donor_disease.to_csv(path_results+'Donor_Disease'+cl_str+'.csv')
        df_donors.to_csv(path_results+'Number_cells_per_donor_and_celltype_'+opt_celltype_groups+cl_str+'.csv')
        df_p_donors.to_csv(path_results+'Percentage_cells_per_donor_and_celltype_'+opt_celltype_groups+cl_str+'.csv')
    elif opt_df_cells_per_donor_and_CT == 'load':
        os.chdir(path_data)
        df_donor_disease = pd.read_csv(path_results+'Donor_Disease'+cl_str+'.csv')
        df_donors = pd.read_csv(path_results+'Number_cells_per_donor_and_celltype_'+opt_celltype_groups+cl_str+'.csv')
        df_p_donors = pd.read_csv(path_results+'Percentage_cells_per_donor_and_celltype_'+opt_celltype_groups+cl_str+'.csv')
        
    if opt_plot:
        ut.plot_cells_per_donor_for_each_cell_type_cluster(df_donors,df_donor_disease,pd.DataFrame(),[],True,path_results,'Number',opt_celltype_groups,n_clusters,cl_str)
        ut.plot_cells_per_donor_for_each_cell_type_cluster(df_p_donors,df_donor_disease,pd.DataFrame(),[],True,path_results,'Percentage',opt_celltype_groups,n_clusters,cl_str)
        
    if opt_downsample_cells:
        np.random.seed(123515)
        #get min percentage of cells across CTs while removing at max 10% of samples = 8 samples (plus any sample that had nan number of cells):
        CTs_assigned = [col for col in df_p_donors.columns if col not in ['Donor', 'Disease']]
        #number of cells sampled for each CT
        #for each CT same number of cells picked
        cutoff_number_of_cells = int(np.floor(np.nanquantile(df_donors[CTs_assigned],0.2,axis=0).min()))
        print('cutoff_number_of_cells: '+str(cutoff_number_of_cells))
        #for a specific cell type: collapse/ aggregate counts from at least 0.2% of cells for each sample (randomly picked), exclude samples that have less than 0.2% of cells
        if 'Unnamed: 0' in CTs_assigned: CTs_assigned.remove('Unnamed: 0')
        n_donors = len(df_p_donors['Donor'].tolist())
        bool_first = True
        for c_id,ct in enumerate(CTs_assigned):
            #for every cell type a different set of donors gets excluded from the analysis
            donors_to_exclude = df_donors[df_donors[ct]<cutoff_number_of_cells]['Donor'].tolist()
            #n_cells_to_pick = int(np.floor(df_donors['Excitatory_10'].sum()/(n_donors-len(donors_to_exclude))*0.2))
            #n_cells_to_pick = int(df_donors[~df_donors['Donor'].isin(donors_to_exclude)][ct].min())
            #store number of cells picked and which donors were excluded in data frame
            if bool_first:
                df_downsampling_summary = pd.DataFrame(data={'Cell type cluster':ct, 'n_cells_selected':cutoff_number_of_cells, 'n_donors_to_exclude':len(donors_to_exclude), 'donors_to_exclude':','.join(donors_to_exclude)},index=[0])
                bool_first = False
            else:
                df_downsampling_summary_add = pd.DataFrame(data={'Cell type cluster':ct, 'n_cells_selected':cutoff_number_of_cells, 'n_donors_to_exclude':len(donors_to_exclude), 'donors_to_exclude':','.join(donors_to_exclude)},index=[0])
                df_downsampling_summary=df_downsampling_summary.append(df_downsampling_summary_add,ignore_index=True,sort=False)
                
            # for CT with fewest cells: how many cells are left per donor?
            # rather take same number of cells from each cell type instead of fixed percentage
            ct_loom_filename = 'Samples_cellranger'+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups+'_'+ct+'.loom'
            for i in range(0,n_datasets):
                ct_downsampled_loom_filename = 'Samples_cellranger'+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups+'_'+ct+'_downsampled_'+str(i)+'.loom'
                # randomly pick n_cells_to_pick from each donor
                os.chdir(path_data_source)
                with loompy.connect(ct_loom_filename) as D_CT:
                    Donor_list = np.unique(D_CT.ca['Donor'])
                    for d_ex in donors_to_exclude:
                        if d_ex in Donor_list:
                            Donor_list=np.delete(Donor_list,np.where(Donor_list==d_ex)[0][0])
                    if len(Donor_list)>40:
                        cells_kept = np.zeros((1, np.shape(D_CT)[1]), dtype=bool).flatten()
                        for donor in Donor_list:
                            IDs = np.where(D_CT.ca['Donor']==donor)
                            IDs_sel = np.random.choice(np.asarray(IDs)[0],cutoff_number_of_cells, replace=False)
                            cells_kept[np.sort(IDs_sel)]=np.ones((1, len(IDs_sel)), dtype=bool).flatten()
                        genes_kept = np.ones((1, np.shape(D_CT)[0]), dtype=bool).flatten()
                        if np.sum(cells_kept)>0:
                            ut.create_loom_file(D_CT, genes_kept, cells_kept, path_data, ct_downsampled_loom_filename)
        if opt_plot:
            ut.plot_cells_per_donor_for_each_cell_type_cluster(df_donors,df_donor_disease,df_downsampling_summary,[],True,path_results,'Number',opt_celltype_groups,n_clusters,cl_str)
            ut.plot_cells_per_donor_for_each_cell_type_cluster(df_p_donors,df_donor_disease,pd.DataFrame(),cutoff_number_of_cells,True,path_results,'Percentage',opt_celltype_groups,n_clusters,cl_str)
   