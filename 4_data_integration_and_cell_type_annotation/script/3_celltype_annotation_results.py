# -*- coding: utf-8 -*-
#date: Created on Fri Nov 27 17:19:57 2020
#@author: Lisa Bast
# version: 0.2
#about: integrate and visualize results from scmap, UMAP with metadata Fig. S1 G,H

import os
import numpy as np
import pandas as pd
import loompy
import sys
import glob
#import seaborn as sns
#from matplotlib import cm

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

main_path = path_code.replace("4_data_integration_and_cell_type_annotation\\script","")

opt_server=True
opt_D_version=''#'v1'#'v2'
ref_data_name = 'Allen_Brain_Human_Multiple_Cortical_Areas_SMARTseq' #'Allen_Brain_Human_M1_10X' 
ref_data_name_short = 'ABM_MCA' #'ABM_M1'
opt_n_CT = 51#76

opt_ref_data_as_test_data=False#True#
opt_load_annotaion_DF = False#True#
opt_plot = True#False#

if opt_plot:
    opt_clustering_plots=False#True#
    opt_plot_annotation_percentages_and_scores = False#True#
    opt_clustering_plots_single_sample = False#True #
    opt_create_cellclass_specific_files = False#True
    opt_clustering_plots_of_cell_classes=True#False#
    opt_plot_cluster_purity_results = True
    opt_abundance = "relative"
    opt_sep_4_scz_ctrl = False # separate plots for the groups (SCZ, CTRL)
    opt_plot_umap = True
else:
    opt_clustering_plots=False
    opt_plot_annotation_percentages_and_scores = False
    opt_clustering_plots_single_sample = False
    opt_clustering_plots_of_cell_classes=False
    opt_plot_cluster_purity_results = False

if opt_clustering_plots_of_cell_classes or opt_create_cellclass_specific_files:
    opt_use_subsampled_data=False
    if opt_use_subsampled_data:
        target_percentage_of_cells=20
    else:
        target_percentage_of_cells=[]
    #opt_cell_classes='scmap'
    opt_cell_classes='conos'

n_genes_ref = [500,1000,2000,3000,5000]
red_methods = ['HDO','HVG_with_seurat_v3','HVG_with_cell_ranger']

if opt_ref_data_as_test_data==True:
    add_str = '_test'
    add_str2 = 'Test_'
else:
    add_str = ''
    add_str2 = 'Filtered_'
    #opt_integrate_annotation_in_loom_file_metadata = False#True#
    opt_merge_annotated_loom_files = True
    opt_integrate_conos_cluster_in_loom_file_metadata = False#True#
    if opt_integrate_conos_cluster_in_loom_file_metadata:
        opt_conos_resolution="7"
        opt_merge_samples=True

tool = 'scmap'
if tool =='scmap' and (opt_plot_annotation_percentages_and_scores):
    #settings=['cluster','2000','HVG_with_cell_ranger_']
    settings=['cell2cluster','5000','HVG_with_cell_ranger_']
    
add_str_align= '_cellranger'

path_query_data, _ , path_results = ut.get_paths(main_path,"Cell type annotation")
path_results = path_results+'/'+tool+'/'+add_str2+ref_data_name_short+'/'+add_str_align[1:]+'/'+str(opt_n_CT)+'_CTs/'
path_ref_data = main_path + '/4_data_integration_and_cell_type_annotation/output/reference_data_sets/'+ref_data_name_short+'/'

path_results_figures = path_results + 'figures_'+opt_D_version+'/'
path_results_files = path_results + 'output_files_'+opt_D_version+'/'
cwd = os.getcwd()

if opt_ref_data_as_test_data:
    # ground truth:
    os.chdir(path_ref_data)

    GT = pd.read_csv('matrix_2500_cells_GT_test.csv')
    GT.rename(columns = {'Unnamed: 0':'Original_Cell_ID'},inplace=True)
    #add cell class
    GT = ut.add_cell_class_info(GT,'cell type')
    
    #create dataframe of results
    os.chdir('./red/')
    ann_all = pd.read_csv('9606_map.csv')
    ann_all.rename(columns = {'Unnamed: 0':'Original_Cell_ID'},inplace=True)
    DF_R = ut.get_DF_annotation_result_test(GT, red_methods, tool, n_genes_ref, path_results_files)
    
    #plot annotation success
    if opt_plot:
        variables = ['P_correct','P_unassigned','P_false','Cohens_kappa','P_CC_correct','P_CC_unassigned','P_CC_false']
        ut.plot_annotation_success_percentage(DF_R,variables,path_results_figures)

    #print best performing methods
    ut.print_best_performing_CT_annotation_modes(DF_R)
    #print difference between best performing methods
    ut.print_difference_between_best_performing_methods(DF_R)
    
    ut.plot_UMAP_test_train('matrix_split_into_test_and_train.loom', path_ref_data, path_results_figures)

else:
    if opt_merge_annotated_loom_files:
        target_files = glob.glob('S*'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered.loom')
        target_files_new = [t for t in target_files if not t.startswith('Samples')]
        loompy.combine(files = target_files_new, output_file = 'Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom', batch_size = 3000)
    if opt_integrate_conos_cluster_in_loom_file_metadata:
        path_loom_file = main_path+"data/filtered_loom_formatted_data"+add_str_align+"/"
        loom_file_name = 'Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom'
        os.chdir(main_path+"data/filtered_loom_formatted_data"+add_str_align+"/")
        df_c = pd.read_csv("conos_clusters_res_"+opt_conos_resolution+".csv",sep=";",names=['Cell_ID','conos_cluster'])
        samples_to_skip = []
        ut.integrate_conos_cluster_in_metadata(df_c,path_loom_file,loom_file_name, opt_merge_samples,opt_conos_resolution,samples_to_skip)
    if opt_plot:
        if opt_plot_annotation_percentages_and_scores:
            os.chdir(path_query_data)
            ca_ann = ['CT_ann_'+ref_data_name_short+'_'+tool+'_'+settings[0]+'_'+str(opt_n_CT)+'CTs']
            ca_score = ['CT_ann_score_'+ref_data_name_short+'_'+tool+'_'+settings[0]+'_'+str(opt_n_CT)+'CTs']
            D = loompy.connect('Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom')
            #ca_ann = ['CT_ann_ABM_MCA_scmap_cluster_'+str(opt_n_CT)+'CTs']#['CT_ann_ABM_MCA_scmap_cell2cluster', 'CT_ann_ABM_MCA_scmap_cluster']
            #ca_score = ['CT_ann_scoreABM_MCA_scmap_cluster_'+str(opt_n_CT)+'CTs']
            sample_IDs,n_cells,disease_status,sex_status = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D)
            sample_IDs_CTRL = [i for i,d in enumerate(disease_status) if d=='CTRL']
            sample_IDs_SCZ = [i for i,d in enumerate(disease_status) if d=='SCZ']
            if opt_load_annotaion_DF==False:
                DF_CT_A, UA_per_donor = ut.get_DF_annotation_result(D,ca_ann,ca_score,sample_IDs,n_cells,disease_status,sex_status)
                DF_CT_A.to_pickle(path_results_files+"DF.pkl")
                np.save(path_results_files+'Unassigned_cells_per_donor.npy',UA_per_donor)
            else:
                DF_CT_A = pd.read_pickle(path_results_files+"DF.pkl")
                UA_per_donor = np.load(path_results_files+'Unassigned_cells_per_donor.npy')
            
            
            #1) plot % of cells unassigned per sample (barplot, color of bar referring to scz/ ctrl)  
            ut.plot_percentage_unassigned_per_sample_and_method(UA_per_donor,ca_ann,sample_IDs,sample_IDs_CTRL,sample_IDs_SCZ,path_results_figures)
            #2) plot % of cells assigned to cell type x ((a) x subplots with barplots having samples on the x axis, color coded by scz/ ctrl),
            #                                            (b) heatmap samples vs cell types with color gradient for range (0%,max%)
            #                                            (c) PCA on cell type composition
            ut.plot_percentage_assigned_to_each_CT_per_method(ca_ann,DF_CT_A,path_results_figures)
            ut.plot_heatmap_CT_vs_sample(DF_CT_A,ca_ann,path_results_figures, 'percentage_assigned_to_each_CT_per_sample')
            ut.plot_PCA_CT_percentages(DF_CT_A,ca_ann,path_results_figures)
            
            #3) heatmap samples vs cell types with mean/ median score as color label
            ut.plot_heatmap_CT_vs_sample(DF_CT_A,ca_ann,path_results_figures, 'median_score_assigned_to_each_CT_per_sample')
            ut.plot_heatmap_CT_vs_sample(DF_CT_A,ca_ann,path_results_figures, 'mean_score_assigned_to_each_CT_per_sample')
            #4) cluster purity: for each cluster (leiden) calculate % of most abundant cell type
    
        if opt_clustering_plots==True:
            #if opt_plot_UMAP_of_cell_classes==False:
            #UMAP with assigned cell types as color label
            #UMAP with score as color label
            single_marker_gene_list,single_marker_CT_list = ut.get_marker_list()
            if opt_clustering_plots_single_sample:
                add_str_figure_file_name = 'S2'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated'
                leiden_res = [0.5,0.6]
                conos_res=[7]
                #variables_color_code = ['CT_ann_ABM_MCA_scmap_cell2cluster','CT_ann_ABM_MCA_scmap_cluster','CT_ann_score_ABM_MCA_scmap_cell2cluster','CT_ann_score_ABM_MCA_scmap_cluster','Disease', 'Donor', 'Library', 'PMI_h', 'Sex', 'Age', 'mean_counts_per_barcode', 'mean_reads_per_umi', 'median_gpc', 'num_reads', 'num_umi', 'p_cell_del_filt', 'p_genome_not_gene', 'p_mapped_reads', 'p_sequencing_saturation', 'p_unmapped_reads', 'p_valid_barcodes', 'std_counts_per_barcode']
                variables_color_code = ['CT_ann_ABM_MCA_scmap_cluster_76CTs','CT_ann_scoreABM_MCA_scmap_cluster_76CTs','CT_ann_ABM_MCA_scmap_cluster_51_CTs','CT_ann_score_ABM_MCA_scmap_cluster_51_CTs','Conos_cluster_res_7']
            else:
                add_str_figure_file_name = 'Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated'
                leiden_res = [1.5]
                conos_res=[7]
                if opt_server == True:
                    variables_color_code = ['Conos_cluster_res_7', 'CT_ann_ABM_MCA_scmap_cell2cluster_76CTs', 'CT_ann_score_ABM_MCA_scmap_cell2cluster_76CTs', 'Disease', 'Donor', 'Library', 'Sex', 'Age','PMI_h', 'mean_counts_per_barcode', 'mean_reads_per_umi', 'median_gpc', 'num_reads', 'num_umi', 'p_cell_del_filt', 'p_genome_not_gene', 'p_mapped_reads', 'p_sequencing_saturation', 'p_unmapped_reads', 'p_valid_barcodes', 'std_counts_per_barcode']
                else:
                    #variables_color_code = ['CT_ann_ABM_MCA_scmap_cell2cluster','CT_ann_ABM_MCA_scmap_cluster','CT_ann_score_ABM_MCA_scmap_cell2cluster','CT_ann_score_ABM_MCA_scmap_cluster','Disease', 'Donor', 'Library', 'PMI_h', 'Sex', 'Age', 'mean_counts_per_barcode', 'mean_reads_per_umi', 'median_gpc', 'num_reads', 'num_umi', 'p_cell_del_filt', 'p_genome_not_gene', 'p_mapped_reads', 'p_sequencing_saturation', 'p_unmapped_reads', 'p_valid_barcodes', 'std_counts_per_barcode']
                    variables_color_code = ['Conos_cluster_res_7', 'CT_ann_ABM_MCA_scmap_cell2cluster_76CTs','CT_ann_score_ABM_MCA_scmap_cell2cluster_76CTs','Disease', 'Donor', 'Library', 'Sex']#, 'Age','PMI_h', 'mean_counts_per_barcode', 'mean_reads_per_umi', 'median_gpc', 'num_reads', 'num_umi', 'p_cell_del_filt', 'p_genome_not_gene', 'p_mapped_reads', 'p_sequencing_saturation', 'p_unmapped_reads', 'p_valid_barcodes', 'std_counts_per_barcode']
            loom_file_name = add_str_figure_file_name+'.loom'
            opt_which_clustering='conos_and_leiden'
            ut.plot_PCA_UMAP_clusterContribution_from_loom(loom_file_name, path_query_data, add_str_figure_file_name,path_results,variables_color_code,single_marker_gene_list,single_marker_CT_list,leiden_res,opt_server,opt_clustering_plots_single_sample,opt_which_clustering,[],opt_abundance,opt_sep_4_scz_ctrl,opt_plot_umap)

            if opt_clustering_plots_of_cell_classes:
                if 'D' in locals() or 'D' in globals():
                    D.close()
                #UMAP plots for 3 classes inhibitory/gabaergic, exhibitory/glutamatergic, non_neuronal
                CT_classes = ['Excitatory','Inhibitory','NonNeuronal']
                if opt_create_cellclass_specific_files==True:
                    file_names_new = ut.create_cellclass_specific_files(path_query_data,CT_classes,opt_D_version,opt_cell_classes,opt_use_subsampled_data,target_percentage_of_cells,opt_n_CT)
                else:
                    file_names_new = ['Samples_Excitatory_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom',
                                      'Samples_Inhibitory_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom',
                                      'Samples_NonNeuronal_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom']
                os.chdir(path_query_data)
                #with loompy.connect('Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated.loom') as D:
                for f_id,fnn in enumerate(file_names_new):
                    opt_which_clustering = 'conos_and_leiden'
                    leiden_res = [0.5,0.6]
                    conos_res=[7]
                    opt_clustering_plots_single_sample=False
                    if opt_plot and opt_clustering_plots_of_cell_classes:
                        ut.plot_PCA_UMAP_clusterContribution_from_loom(fnn, path_query_data, fnn[:-5], path_results_figures,variables_color_code,single_marker_gene_list,single_marker_CT_list,leiden_res,opt_server,opt_clustering_plots_single_sample,opt_which_clustering,[],opt_abundance,opt_sep_4_scz_ctrl,opt_plot_umap)

        if opt_plot_cluster_purity_results:
            #visualize cluster purity results
            os.chdir(path_results)
            leiden_res = 1.5
            conos_res=7
            if opt_clustering_plots_of_cell_classes:
                folder_names = ['Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated',
                                'Samples_Excitatory_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated',
                                'Samples_Inhibitory_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated',
                                'Samples_NonNeuronal_'+opt_cell_classes+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated']
            else: 
                folder_names = ['Samples'+add_str_align+'_TH_and_D_adj_'+opt_D_version+'_filtered_and_CT_annotated']
            file_names = ['P_cells_CT_ann_ABM_MCA_scmap_cluster_'+str(opt_n_CT)+'CTs_leiden_'+str(leiden_res)+'_.xlsx'] #'P_cells_CT_ann_ABM_MCA_scmap_cluster_'+str(opt_n_CT)+'CTs_Conos_cluster_res_'+str(conos_res)+'_.xlsx']#,
            clustering_mode = 'conos_res_'+str(conos_res)#['leiden_res_'+str(leiden_res)]#
            for f_id,folder_name in enumerate(folder_names):
                os.chdir(path_results+folder_name)
                for file_id,file_name in enumerate(file_names):
                    df = pd.read_excel(file_name, sheet_name='raw') 
                    df_max_per_cluster = pd.read_excel(file_name, sheet_name='max_per_cluster')
                    df_max_per_celltype = pd.read_excel(file_name, sheet_name='max_per_celltype')
                    ut.plot_hist_max_contribution_per_cluster(df_max_per_cluster,path_results+folder_name+'/figures_'+opt_D_version+'/',clustering_mode[file_id])
                    #fp.plot_CT_specificity_per_cluster
                
                    df_max_per_celltype.rename(columns={0:'max contribution to a cluster'},inplace=True)
                    df_max_per_celltype['max contribution to a cluster'] = df_max_per_celltype['max contribution to a cluster']*100
                    df_max_per_cluster.drop(columns=0,inplace=True)

#plot heatmap with seurat using loom files in R
        

        