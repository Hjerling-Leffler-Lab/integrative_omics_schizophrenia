# about:  plot results of DEG (Differentially expressed genes) analysis scRNAseq/ snRNAseq data 
#         performed with DESeq2 
#         Fig. 2 B-E,H, S6 A-E
# Created on Wed Mar  3 14:44:47 2021
# @author: Lisa Bast
# version: 0.3.2

import anndata
import scanpy as sc
import os
import numpy as np
import pandas as pd
import sys
import statsmodels.stats.multitest as ssm
import gc
    
## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("5b_differential_gene_expression_analysis\\script","")

if __name__ == "__main__": ut.clear_ws()
gc.collect()

#specify settings

#define colors
col_up_DEG='darkred'
col_down_DEG = 'teal'

opt_pagoda_filtering=True
opt_save=True
#opt_celltype_groups = 'scmap'
opt_celltype_groups = 'conos_cluster_based'
opt_plot_number_deregulated_genes = True#False#
opt_plot_correlation_clustermap = False#True
opt_plot_volcano = False#True#
opt_plot_log2FC_per_CT = True#False#
opt_comparison_with_genelists = False#True#
opt_add_proteomics_to_genelist_comparison = True
opt_analyze_correlation_with_genelists=False#True#False#
alpha_val_range=[0.3,0.05]
opt_export_genes_ensemble = True#False #
cutoff_q_p=0.3#gene lists in genematrix

if opt_comparison_with_genelists or opt_analyze_correlation_with_genelists:
    genelist_version = 'v2' #'v1'#'mildas_panel' #
    opt_genes_to_show = 'top_genes'# 'genes_with_strong_signal' #
    if opt_genes_to_show=='top_genes':
        top_n_genes=50
        opt_par = top_n_genes
    elif opt_genes_to_show =='genes_with_strong_signal':
        cutoff = 0.2
        opt_par = cutoff

add_str = '_cellranger'
n_cluster=[16,37]#,3]

if opt_pagoda_filtering:
    add_str_pagoda = '_pagoda'

correction_methods=["benjamini-hochberg","bonferroni"]

DESeq2_design_mode = 'with_6_RUVs'
DESeq2_aggregation_methods =['sum']

path_results_proteomics = path_project+"5b_differential_gene_expression_analysis/data/proteomics/" #old: path_project+'data/Gene_lists/'
path_results_proteomics_hgnc = path_project+"5b_differential_gene_expression_analysis/data/proteomics/"
hgnc_file = "hgnc_2022-04-24_small.tsv"
proteomics_file = "diagnosis_vs_layers.xlsx" #differential_abundance_analysis.xlsx"

path_filtered_loom_data = path_project + '4_data_integration_and_cell_type_annotation/output/'
gene_matrix_file = path_project + "5b_differential_gene_expression_analysis/data/Gene_lists/geneMatrix_v2.tsv"
bool_first_alpha = True
for alpha_val in alpha_val_range:
    for n_cl in n_cluster:
        # read data as anndata object
        #investigate results
        path_data, path_results = ut.get_paths(path_project,'DEG_visualization')
        if n_cl!=3:
            path_results_short = path_results + str(n_cl) + "_CTs"
        else:
            path_results_short = path_results + str(n_cl) + "_classes"
        #create list of cell types:
        os.chdir(path_data+'Aggregated_data/')
        str_start = 'Aggregated_counts'+add_str+add_str_pagoda+'_TH_and_D_adj_filtered_'+opt_celltype_groups
        csv_aggr_files = [ f for f in os.listdir( os.curdir ) if os.path.isfile(f) and f.startswith(str_start) and f.endswith('.csv')]
        CTs_assigned = [f.split("_based_")[1].split('_'+DESeq2_aggregation_methods[0]+'.csv')[0] for f in csv_aggr_files]
        
        genes_excitatroy=[];
        path_results_long = path_results_short + 'design_'+DESeq2_design_mode+'/'
        os.chdir(path_results_long)
        for j,am in enumerate(DESeq2_aggregation_methods):
            celltype_folders = [ f.path for f in os.scandir(path_results_long) if f.is_dir()]
            bool_first_DF_DESeq2 = True
            bool_first_pvals = True
            bool_first_CT = True
            bool_first_df_true = True
            bool_first_df_rand_sign = True
            method_str = DESeq2_design_mode+'_'+am
            for k,ctf in enumerate(celltype_folders):
                os.chdir(ctf)
                #get adjusted p-values for all genes tested:
                if os.path.isfile('df_results_no_shrinkage.csv'):
                    df_DESeq2_all_genes = pd.read_csv('df_results_no_shrinkage.csv',sep=',',names=["Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"])
                    df_DESeq2_all_genes.drop(index=0, inplace=True)
                    df_DESeq2_all_genes = df_DESeq2_all_genes.astype({"baseMean":'float',"log2FoldChange":'float',"lfcSE":'float',"stat":'float',"pvalue":'float',"padj":'float'})
                else:
                    df_DESeq2_all_genes = df_DESeq2_all_genes.rename(columns={"Unnamed: 0": "Gene"})
                    continue
                CT_i = "_".join(ctf.split('/')[len(ctf.split('/'))-1].split('_')[1:-1])
                p_vals_i = df_DESeq2_all_genes['padj']
                ##TO DO: hier weiter: for each CT add if gene up or down
                bool_up_regulated = df_DESeq2_all_genes['log2FoldChange'] > 0
                #rename column of genes:
                
                #split in gene and ensgid (easier for GSA):
                genes = df_DESeq2_all_genes['Gene'].tolist()#df_DESeq2_all_genes["Gene"]
                genes_short = [s.split('_')[0] for s in genes]#[s.split('_')[0] for s in df_DESeq2_all_genes["Gene"].tolist()]
                ensgids = [s.split('_')[1] for s in genes]#[s.split('_')[1] for s in df_DESeq2_all_genes["Gene"].tolist()]
                #add genes_short as column and make index:
                df_DESeq2_all_genes['Gene_short'] = genes_short
                df_DESeq2_all_genes.set_index('Gene_short',inplace=True)
                #add log10 of adj p-vals as column:
                df_DESeq2_all_genes['negative_log10_padj'] = -np.log10(df_DESeq2_all_genes['padj'])
                #add color for sign genes as column:
                df_DESeq2_all_genes['Color_sign'] = ['grey']*np.shape(df_DESeq2_all_genes)[0]
                df_DESeq2_all_genes.loc[(df_DESeq2_all_genes.padj<=alpha_val) & (df_DESeq2_all_genes.log2FoldChange<=0),'Color_sign']=col_down_DEG
                df_DESeq2_all_genes.loc[(df_DESeq2_all_genes.padj<=alpha_val) & (df_DESeq2_all_genes.log2FoldChange>0),'Color_sign']=col_up_DEG

                df_DESeq2_all_genes_CT_i = df_DESeq2_all_genes.copy()
                df_DESeq2_all_genes_CT_i['celltype'] = CT_i

                #vulcano plot:
                if opt_plot_volcano:
                    ut.plot_volcano(df_DESeq2_all_genes_CT_i,CT_i,path_results_long,alpha_val,opt_save,annotations="very_small_pvals")
                
                if opt_analyze_correlation_with_genelists or opt_export_genes_ensemble or opt_comparison_with_genelists or opt_plot_log2FC_per_CT:
                    if bool_first_CT:
                        df_DESeq2_all_genes_all_CTs = df_DESeq2_all_genes_CT_i.copy()
                        if opt_analyze_correlation_with_genelists:
                            df_corr_genelists_up,df_corr_genelists_down,bool_first_CT = ut.get_comparison_DEG_pvals_to_genelists(df_DESeq2_all_genes,opt_save,path_results_long,CT_i,alpha_val,cutoff_q_p,bool_first_CT,[],[],path_project,genelist_version)
                        bool_first_CT=False
                    else:
                        df_DESeq2_all_genes_all_CTs = pd.concat([df_DESeq2_all_genes_all_CTs,df_DESeq2_all_genes_CT_i])
                        #compare to gene lists from gene matrix:
                        if opt_analyze_correlation_with_genelists:
                                df_corr_genelists_up,df_corr_genelists_down,bool_first_CT = ut.get_comparison_DEG_pvals_to_genelists(df_DESeq2_all_genes,opt_save,path_results_long,CT_i,alpha_val,cutoff_q_p,bool_first_CT,df_corr_genelists_up,df_corr_genelists_down,path_project,genelist_version)
                        
                #make sure to be back to correct directory:
                os.chdir(ctf)
                
                if bool_first_alpha:
                    df_p_vals_i = pd.DataFrame(data={'Gene':genes,'Gene_short':genes_short,'ensgid':ensgids, CT_i : p_vals_i, 'bool_'+CT_i+'_upregulated': bool_up_regulated})
                    if bool_first_pvals:
                        DF_p_vals = df_p_vals_i.copy()
                        bool_first_pvals = False
                    else:
                        DF_p_vals = pd.merge(DF_p_vals,df_p_vals_i,on=["Gene","Gene_short","ensgid"],how='outer')
                        #now contains some nans
                else:
                    DF_p_vals = pd.read_csv(path_results_long+'df_p_vals_adj_'+am+'_per_celltype.csv')
                
                #get results for significant genes:
                if os.path.isfile('df_results_'+am+'.csv'):
                    df_DESeq2 = pd.read_csv('df_results_'+am+'.csv')
                elif os.path.isfile('df_results_'+am+'_RUV.csv'):
                    df_DESeq2 = pd.read_csv('df_results_'+am+'_RUV.csv')
                elif os.path.isfile('df_results_'+am+'_MARS.csv'):
                    df_DESeq2 = pd.read_csv('df_results_'+am+'_MARS.csv')
                else:
                    continue
                df_DESeq2 = ut.add_gene_names(df_DESeq2)
                
                # export all excitatory genes:
                if "Excitatory" in ctf:
                    genes_excitatroy.append(df_DESeq2['genes'].tolist())
                    
                if bool_first_DF_DESeq2:
                    DF_DESeq2 = df_DESeq2.copy()
                    bool_first_DF_DESeq2 = False
                else:
                    DF_DESeq2 = pd.concat([DF_DESeq2, df_DESeq2])

            # last CT
            df_DESeq2_all_genes_all_CTs.to_csv(path_results_long+'T5_DEGs_per_cell_type.csv')

            #plot signed p-value pairwise correlation between cell types as clustermap
            if opt_plot_correlation_clustermap:
                ut.plot_signed_DEG_corr_clustermaps(df_DESeq2_all_genes_all_CTs,path_results_long)

            if alpha_val==0.05:
                #effect size per CT:
                if opt_plot_log2FC_per_CT:
                    ut.plot_log2FC_per_CT(df_DESeq2_all_genes_all_CTs,alpha_val,n_cl,path_filtered_loom_data,path_results_long,opt_save)
            
            #number deregulated genes per cell type:
            if opt_plot_number_deregulated_genes:
                #TO DO: make sure all CTs are in DF
                ut.plot_number_of_deregulated_genes(df_DESeq2_all_genes_all_CTs,'no_shrinkage',opt_save,path_results_long,'','vs_celltypes',n_cl,path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_loom_data,gene_matrix_file,n_cl,True)
                ut.plot_number_of_deregulated_genes(df_DESeq2_all_genes_all_CTs,'no_shrinkage',opt_save,path_results_long,'','vs_number_cells',n_cl,path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_loom_data,gene_matrix_file,n_cl,False)

            if alpha_val==0.05:
                #only interested in this for alpha=0.05
                #add cell type names to DF_DESeq2:
                #DF_DESeq2 = fg.add_celltype_names_to_DF(DF_DESeq2, 'celltype',opt_velo, n_cl, path_filtered_loom_data)
                if opt_analyze_correlation_with_genelists:
                    #save df_corr_genelists
                    df_corr_genelists_up.to_csv(path_results_long+'correlation_of_upreg_DEGs_adj_pval_with_genematrix_lists_qvalues.csv')
                    df_corr_genelists_down.to_csv(path_results_long+'correlation_of_downreg_DEGs_adj_pval_with_genematrix_lists_qvalues.csv')
                #comparison with geneMatrix (venn plot, all DEGs together):
                #load genelists from geneMatrix for comparison 
                if opt_comparison_with_genelists:
                    genelists_wp, genelist_names_wp = ut.get_gene_lists(path_project,genelist_version)
                    if genelist_version == 'v1' or genelist_version == 'v2':
                        for i,gl in enumerate(genelists_wp):
                            ut.plot_venn_geneMatrix(DF_DESeq2,gl,genelist_names_wp[i],opt_save,path_results_long,am,'no_shrinkage','',alpha_val)
            
                    #genelists plots (independent of alpha_val --> only run once)
                    genelists, genelist_names, q_values = ut.get_gene_lists_with_qvalues(path_project,genelist_version)
                    id_last_gl = len(genelist_names)
                    genelists.extend(genelists_wp)
                    genelist_names.extend(genelist_names_wp)
                    n_genes_total = np.zeros(np.shape(genelists))
                    n_genes_plotted = np.zeros(np.shape(genelists))
                    if opt_comparison_with_genelists:
                        for i_gl,gl in enumerate(genelists):
                            if i_gl < id_last_gl:
                                q_vals = q_values[i_gl]
                            else:
                                q_vals=[]
                            n_genes_total[i_gl], n_genes_plotted[i_gl] = ut.plot_heatmap_DEG_result_for_genes_in_genelist(df_DESeq2_all_genes_all_CTs,gl,genelist_names[i_gl],q_vals,opt_save,path_results_long,opt_genes_to_show, opt_par,genelist_version, 'p-value',opt_add_proteomics_to_genelist_comparison,path_results_proteomics, path_results_proteomics_hgnc, hgnc_file, proteomics_file)
                            ut.plot_heatmap_DEG_result_for_genes_in_genelist(df_DESeq2_all_genes_all_CTs,gl,genelist_names[i_gl],q_vals,opt_save,path_results_long,opt_genes_to_show, opt_par,genelist_version, 'log2FC',opt_add_proteomics_to_genelist_comparison,path_results_proteomics, path_results_proteomics_hgnc, hgnc_file, proteomics_file)
                        qval_cutoff_range = [0.05,0.1,0.2,0.3]#,0.5,1]
                        if genelist_version == 'v1' or genelist_version == 'v2':
                            ut.plot_percentage_genes_in_list_and_DEG_for_range_of_cutoffs(path_project, df_DESeq2_all_genes_all_CTs, opt_par, qval_cutoff_range, path_results_long, opt_save)
        
            #for a range of alpha values:
            if opt_export_genes_ensemble:
                #filter for pval below alpha
                ut.get_ensemble_genes(df_DESeq2_all_genes_all_CTs,alpha_val,path_results_long,False)
                if n_cl == 3 or n_cl==15:
                    #per cell type
                    ut.get_ensemble_genes(df_DESeq2_all_genes_all_CTs,alpha_val,path_results_long,True)

            if bool_first_alpha:
                DF_p_vals.to_csv(path_results_long+'df_p_vals_adj_'+am+'_per_celltype.csv')
                bool_first_alpha = False

            #no fdr correction for different alpha_vals
            DF_DESeq2_not_fdr_corr = ut.get_results_DF_for_sign_genes(DF_p_vals,celltype_folders,alpha_val)
            # DF_DESeq2_not_fdr_corr = fg.add_celltype_names_to_DF(DF_DESeq2_not_fdr_corr, 'Cell type',opt_velo, n_cl, path_filtered_loom_data)
            DF_DESeq2_not_fdr_corr.rename(columns={'Cell type':'celltype_names'},inplace=True)
            if alpha_val==0.05:
                ut.plot_number_of_deregulated_genes(DF_DESeq2_not_fdr_corr,'no_shrinkage',opt_save,path_results_long,'','vs_celltypes',n_cl,path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_loom_data,gene_matrix_file,n_cl,True)
            else:
                ut.plot_number_of_deregulated_genes(DF_DESeq2_not_fdr_corr,'no_shrinkage',opt_save,path_results_long,'','vs_celltypes',n_cl,path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_loom_data,gene_matrix_file,n_cl,False)
            ut.plot_number_of_deregulated_genes(DF_DESeq2_not_fdr_corr,'no_shrinkage',opt_save,path_results_long,'','vs_number_cells',n_cl,path_data,add_str,add_str_pagoda,opt_celltype_groups,alpha_val,path_filtered_loom_data,gene_matrix_file,n_cl,False)
            
            #plot frequency of genes that are detected multiple times for fdr corrected list:
            ut.plot_DEGs_CT_overlap_as_upset(DF_DESeq2_not_fdr_corr,opt_save,path_results_long,am, '',alpha_val)
            
            if alpha_val==0.05:
                #plot heatmap: sort individuals by group; only DEG genes
                DEG_genes_list_all = []
                for c_id,ct in enumerate(CTs_assigned):
                    path = path_results_short+'design_'+DESeq2_design_mode+'/based_'+ct+'_sum/'
                    if os.path.isdir(path):
                        os.chdir(path)
                        if DESeq2_design_mode=='without_covariates':
                            filename = 'df_results_sum.csv'
                        else:
                            filename = 'df_results_sum_RUV.csv'
                        if os.path.isfile(filename):
                            DF_deseq2 = pd.read_csv(filename)
                        else:
                            continue
                        DF_deseq2['names'] = DF_deseq2['Unnamed: 0'].str.split("_",n=1,expand=True)[0]
                        DEG_genes_list = DF_deseq2['names'].tolist()
                        DEG_genes_list_all.extend(DEG_genes_list)
                

    
        

