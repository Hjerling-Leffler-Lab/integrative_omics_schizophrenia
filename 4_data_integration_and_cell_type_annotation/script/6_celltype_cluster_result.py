# -*- coding: utf-8 -*-
# date: Created on Thu Apr 15 17:49:55 2021
# version 0.0.1
# author: Lisa Bast
# about:
#   - integrate 2nd level conos clusters in metadata of loom files
#   - split loom file in 3 main subclasses
#   - optionally subsample data --> required for heatmap plot (R)
#   - merge and remove some 2nd level conos clusters and save final cluster info in loom files
#   - graphics:
#     - plot sankey plot of scmap annotation and final clusters 
#     - plot cluster purity of scmap annoation 
#     - plot scmap annotation specificity for each cluster
#   - integrate cluster names in loom file
#   - Fig. S3 C, S4 A-D
import loompy
import pandas as pd
import os
import numpy as np
import sys
from plotly.offline import plot
import plotly.graph_objects as go

## settings and paths: 
## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

main_path = path_code.replace("4_data_integration_and_cell_type_annotation\\script","")

path_data, _, path_results = ut.get_paths(main_path,"Cell type annotation")

add_str_align = '_cellranger'
opt_server = True

loom_file_path=main_path+'/4_data_integration_and_cell_type_annotation/output/'

opt_plot_umap_scdrs = False#True#
opt_laptop = False

opt_D_version=''#'v1'#
opt_pagoda_filtering=True#False#
opt_abundance = 'total' #'relative'
opt_sep_4_scz_ctrl=False # should groups be plotted separately?
opt_plot_umap = True

samples_to_skip = ['S66','S78']
if opt_pagoda_filtering:
    add_str_pagoda = '_pagoda'
else:
    add_str_pagoda = ''
loom_file_name = 'Samples'+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom'

opt_integrate_conos_cluster_2nd_level = False#True#
opt_subsample_data = False#True#
opt_create_celltype_class_specific_files = False#True#
opt_use_subsampled_data = False#True #
opt_integrate_final_conos_clusters = False#True#

opt_integrate_annotated_CT_reduced = False#True#

opt_merge_CT_class_files=False #True#
opt_plot=True#False#
opt_plot_subsampled_data = False#True
opt_plot_per_class = True #False
opt_plot_sankeys = True 
opt_sankey_of_merged_cluster=True#False#
opt_sankey_without_removed_cells=True#False#
opt_sankey_without_unassigned_cells=True#False#

#opt_cell_classes='scmap'
opt_cell_classes ='conos'
#number of cell types annotated with scmap cluster/ cell2cluster
n_CT_range = [51, 76]#

scmap_mode='cell2cluster' # 'cluster'

if opt_plot_subsampled_data:
    subsample_str = "_subsampled_to_10_percent_cells"
else:
    subsample_str = ""

if (opt_integrate_final_conos_clusters or opt_plot ):
    n_clusters_range = [15,37]
if opt_subsample_data or opt_use_subsampled_data:
    target_percentage_of_cells=10
if opt_use_subsampled_data:
    loom_file_name_subsample = loom_file_name[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom'

os.chdir(loom_file_path)

CT_classes = ['Excitatory','Inhibitory','NonNeuronal']
#integrate 2nd level conos clusters in metadata of loom files
if opt_integrate_conos_cluster_2nd_level:
    #integrate 2nd level conos cluster as metadata for all samples and cells
    for i,c in enumerate(CT_classes):
        csv_file = 'conos_clusters_2ndalignment_'+c+'.csv'
        os.chdir(loom_file_path)
        if i==0:
            DFs = pd.read_csv(csv_file,sep=";",names=['Cell_ID','conos_cluster'])
        else:
            DFs = pd.concat([DFs, pd.read_csv(csv_file,sep=";",names=['Cell_ID','conos_cluster'])])
    ut.integrate_conos_cluster_in_metadata(DFs,loom_file_path,loom_file_name,True,opt_pagoda_filtering,'2nd_level',samples_to_skip)
      
#create CT class specific files            
if opt_create_celltype_class_specific_files:
        #split loom file in 3 main subclasses
        if opt_cell_classes =='conos':
            for c_id,c in enumerate(CT_classes):
                csv_file = 'conos_clusters_2ndalignment_'+c+'.csv'
                os.chdir(loom_file_path)
                DF = pd.read_csv(csv_file,sep=";",names=['Cell_ID','conos_cluster'])
                file_name_new = 'Samples_'+c+'_'+opt_cell_classes+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated'
                # if opt_use_subsampled_data:
                #     loom_file_name_subsample = loom_file_name[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom'
                #     loom_file_name_new = file_names_new[c_id]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom'
                #     D = loompy.connect(loom_file_name_subsample)
                # else:
                loom_file_name_new = file_name_new+'.loom'
                D = loompy.connect(loom_file_name)
                if c_id==0:
                    genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
                A=D.ca['CellID']
                B=DF['Cell_ID'].tolist()
                b=set(B)
                cells_kept = [x in b for x in A]
                ut.create_loom_file(D, genes_kept, cells_kept, loom_file_path,loom_file_name_new)
        elif opt_cell_classes=='scmap':
            os.chdir(loom_file_path)
            if 'D' in locals() or 'D' in globals():
                D.close()
            with loompy.connect(loom_file_name) as D:
                for c_id,c in enumerate(CT_classes):
                    file_name_new = 'Samples_'+c+'_'+opt_cell_classes+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated'+'.loom'
                    if c_id==0:
                        ca_variable = 'CT_ann_ABM_MCA_scmap_'+scmap_mode+'_'+str(n_CT_range[1])+'_CTs'
                        genes_kept = np.ones((1, np.shape(D)[0]), dtype=bool).flatten()
                    cells_kept = np.zeros((1, np.shape(D)[1]), dtype=bool).flatten()
                    metadata_col = ca_variable
                    if c.startswith('Samples_Excitatory'):
                        cells_kept_ids = np.where(np.char.startswith(D.ca[metadata_col].astype(str),'Exc ',start=0,end=None))
                    elif c.startswith('Samples_Inhibitory'):
                        cells_kept_ids = np.where(np.char.startswith(D.ca[metadata_col].astype(str),'Inh ',start=0,end=None))
                    elif c.startswith('Samples_NonNeuronal'):
                        cells_kept_ids = np.where(np.logical_and(np.logical_and(np.logical_and(np.char.startswith(D.ca[metadata_col].astype(str),'Exc ',start=0,end=None)==False,np.char.startswith(D.ca[metadata_col].astype(str),'Inh ',start=0,end=None)==False),np.char.startswith(D.ca[metadata_col].astype(str),'nan',start=0,end=None)==False) , np.char.startswith(D.ca[metadata_col].astype(str),'unassigned',start=0,end=None)==False))
                    cells_kept[cells_kept_ids]=True
                    ut.create_loom_file(D, genes_kept, cells_kept,loom_file_path,file_name_new)
       
#makes changes to ca_variable and stores in new variable
if opt_integrate_annotated_CT_reduced:
    for n_CT in n_CT_range:
        ca_variable = 'CT_ann_ABM_MCA_scmap_'+scmap_mode+'_'+str(n_CT)+'_CTs'
        for CT_class in CT_classes:
            loom_file_name_CT_class = loom_file_name[0:8]+CT_class+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom'
            os.chdir(loom_file_path)
            D = loompy.connect(loom_file_name_CT_class)
            CT_red_col_attr = D.ca[ca_variable]
            inh,exc,nn,unassigned,nans = fg.get_percentages_per_CT_class(D,ca_variable)
            ID_inh,ID_ex,ID_nn = fg.get_IDx_for_CTs_per_CT_class(D,ca_variable)
            if CT_class =='NonNeuronal':
                    CT_red_col_attr[ID_inh]='other'
                    CT_red_col_attr[ID_ex]='other'
            elif CT_class=='Excitatory':
                    CT_red_col_attr[ID_inh]='other'
                    CT_red_col_attr[ID_nn]='other'
            elif CT_class =='Inhibitory':
                    CT_red_col_attr[ID_ex]='other'
                    CT_red_col_attr[ID_nn]='other'
            ut.add_col_metadata_to_loom_file(D, CT_red_col_attr, ca_variable+'_red', loom_file_path, loom_file_name_CT_class, False)
            
#merge and remove some 2nd level conos clusters and save final cluster info in loom files
#integrate final cluster info in cell type class files (D.ca['cluster'],D.ca['cluster_name'])
if opt_integrate_final_conos_clusters:
    for CT_class in CT_classes:
        loom_file_name_CT_class = loom_file_name[0:8]+CT_class+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom'
        # if opt_use_subsampled_data:
        #     loom_file_name_CT_class = loom_file_name_CT_class[0:-5]+'_subsampled_to_'+str(target_percentage_of_cells)+'_percent_cells.loom'
        new_filename = loom_file_name_CT_class
        for n_clusters in n_clusters_range:
            D = loompy.connect(loom_file_name_CT_class)
            new_col_attr, new_col_attr_name = fg.get_final_conos_clustering(D.ca,CT_class,n_clusters)
            opt_automatically_adapt_file_name=False
            if new_col_attr_name != '':
                new_filename = ut.add_col_metadata_to_loom_file(D, new_col_attr, new_col_attr_name, loom_file_path, new_filename, opt_automatically_adapt_file_name)
                #new_filename is overwritten if opt_automatically_adapt_file_name is true
    
                #loom file connection was closed when other loom file was created
                D = loompy.connect(new_filename)
                new_name_col_attr, new_name_col_attr_name, new_color_col_attr, new_color_col_attr_name = fg.get_celltype_names_and_colors_for_clusters(D.ca,n_clusters, loom_file_path)
                new_filename = ut.add_col_metadata_to_loom_file(D, new_name_col_attr, new_name_col_attr_name,loom_file_path, new_filename, opt_automatically_adapt_file_name)
                #new_filename is overwritten if opt_automatically_adapt_file_name is true
                
                D = loompy.connect(new_filename)
                new_filename = ut.add_col_metadata_to_loom_file(D, new_color_col_attr, new_color_col_attr_name, loom_file_path, new_filename, opt_automatically_adapt_file_name)
    
                #CTs_assigned = fg.get_CT_assigned(D.ca,opt_celltype_groups)
            
if opt_merge_CT_class_files:
    loom_file_names_CT_classes=[]
    for CT_class in CT_classes:
        loom_file_names_CT_classes.append(loom_file_name[0:8]+CT_class+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom')
    loompy.combine(files = loom_file_names_CT_classes, output_file = 'Samples'+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom', batch_size = 3000)
    
#optionally subsample data --> required for heatmap plot (R)
if opt_subsample_data:
    # subsample 20% of cells:
    file_names = ['Samples_conos'+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated_and_CT_clustered.loom',
                  loom_file_name]
    for ct_class in CT_classes:
        add_file = 'Samples_'+ct_class+'_'+opt_cell_classes+add_str_align+add_str_pagoda+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated_and_CT_clustered.loom'
        file_names.append(add_file)

    for file in file_names:
       loom_file_name_subsample_class_i = ut.create_subsample_file_percentage_cells(loom_file_path, file, target_percentage_of_cells, False) 
    
if opt_plot:
    #graphics:
     # 1) plot sankey plot of scmap annotation and final clusters
    
    for i,n_CT in enumerate(n_CT_range):
        ca_variable = 'CT_ann_ABM_MCA_scmap_'+scmap_mode+'_'+str(n_CT)+'_CTs'
        ca_variable_score = 'CT_ann_score_ABM_MCA_scmap_'+scmap_mode+'_'+str(n_CT)+'_CTs'
        path_results_figures=main_path+'/4_data_integration_and_cell_type_annotation/output/CT_anno_scmap/ABM_MCA/'
        if opt_plot_sankeys:
            for CT_class in CT_classes:
                lf_name = loom_file_name[0:8]+CT_class+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom'
                print(lf_name)
                os.chdir(loom_file_path)
                D = loompy.connect(lf_name)

                #create a contigency table for cluster and celltype annotation
                for n_clusters in n_clusters_range:
                    #get colors for sankey:
                    if opt_sankey_of_merged_cluster:
                        df_ct_mapping = pd.DataFrame(data={"celltype":D.ca["cluster_name_"+str(n_clusters)+"CTs"],"cluster":D.ca["cluster_"+str(n_clusters)+"CTs"]})
                    else:
                        df_ct_mapping = pd.DataFrame(data={"celltype":D.ca["cluster_name_"+str(n_clusters)+"CTs"],"cluster":D.ca["Conos_2nd_level"]})
                    df_ct_mapping.drop_duplicates(keep="first", inplace=True)
                    df_ct_mapping = ut.add_CT_colors_to_DF(df_ct_mapping,loom_file_path,True,n_clusters)
                    ut.plot_Sankey_of_CT_clusters(CT_class,D,ca_variable,n_clusters,opt_sankey_without_removed_cells, opt_sankey_without_unassigned_cells,opt_sankey_of_merged_cluster,df_ct_mapping,opt_D_version,path_results_figures)
    
        # 2) plot cluster purity of scmap annoation 
        opt_which_clustering = 'final_conos_based'
        single_marker_gene_list,single_marker_CT_list = fg.get_marker_list()
        variables_color_code = ['Disease', 'Donor', 'Library', 'PMI_h', 'Sex', 'Age']
        if i==0:
            variables_color_code =  variables_color_code + [ca_variable,ca_variable_score]
                                    #'mean_counts_per_barcode', 'mean_reads_per_umi', 'median_gpc', 'num_reads', 'num_umi', 'p_cell_del_filt', 'p_genome_not_gene', 'p_mapped_reads', 'p_sequencing_saturation', 'p_unmapped_reads', 'p_valid_barcodes', 'std_counts_per_barcode'
        else:
            variables_color_code =  variables_color_code + [ca_variable] + [ca_variable_score] 
            for n_clusters in n_clusters_range:
                variables_color_code =  variables_color_code + ['cluster_name_'+str(n_clusters)+'CTs'] 
            if opt_plot_per_class:
                for CT_class in CT_classes:
                    lf_name = loom_file_name[0:8]+CT_class+'_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered.loom'
                    add_str_figure_file_name=CT_class+'_'+opt_which_clustering+'_'
                    #conos_cluster_name = 'cluster_'+str(n_clusters)+'CTs'
                    lei_res = [0.4,0.6,0.8]
                    ut.plot_PCA_UMAP_clusterContribution_from_loom(lf_name, loom_file_path, add_str_figure_file_name, path_results_figures,variables_color_code,single_marker_gene_list,single_marker_CT_list,lei_res,opt_server,False,opt_which_clustering,n_clusters_range,opt_abundance,opt_sep_4_scz_ctrl,opt_plot_umap)
            #all classes together:
            print('made it to the final lines of code! Variables color code are:')
            if opt_plot_umap_scdrs:
                lf_name = 'Samples_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered'+subsample_str+'_scdrs.loom'
                variables_color_code =  variables_color_code + ['PASS_Schizophrenia_Pardinas2018'] + ['PASS_Schizophrenia_Singh2022_P_weights'] + ['PASS_Schizophrenia_Singh2022_Zscore_weights']  
            lf_name = 'Samples_'+opt_cell_classes+loom_file_name[7:-5]+'_and_CT_clustered'+subsample_str+'.loom'
            add_str_figure_file_name='Samples_'+opt_which_clustering+'_'
            lei_res = [0.8,0.9,1.0]
            print(variables_color_code)
            ut.plot_PCA_UMAP_clusterContribution_from_loom(lf_name, loom_file_path, add_str_figure_file_name, path_results_figures,variables_color_code,single_marker_gene_list,single_marker_CT_list,lei_res,opt_server,False,opt_which_clustering,n_clusters_range,opt_abundance,opt_sep_4_scz_ctrl,opt_plot_umap)

    # 3) plot scmap annotation specificity for each cluster
    
#how many cells in removed clusters 
#fn ='Samples_Inhibitory_conos_velocyto_TH_and_D_adj_v1_filtered_and_CT_annotated_and_CT_clustered.loom'
#D = loompy.connect(fn)
#len(D.ca['CT_ann_ABM_MCA_scmap_cell2cluster_76CTs'][D.ca['cluster']=='remove'].tolist())
