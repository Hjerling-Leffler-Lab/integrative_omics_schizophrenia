# author: Lisa Bast
# date: 2024-07-04,  13:28:24
# version: 0.0.1
# about: compare DEGs in Ruzika and Batiuk papers with our DEGs, Figure S24

import sys
import numpy as np
import os

## specify project path and load functions:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_project = path_code.replace("11_comparison_with_other_snRNAseq_studies\\script","")
path_main = path_code.replace("script","")
path_results = path_main + "/output/"

opt_load_merged_DFs = False #True#
opt_order = "unsorted" #"cluster_map" #
opt_subset_sign_BAST_genes =True 

file_DEGs_Ruzicka = "D:/Documents/projects/SCZ_human_cortex/data/Gene_lists/Ruzicka_et_al_2024/science.adg5136_data_s4.xlsx" #path_main + "data/Ruzicka_et_al_2024/science.adg5136_data_s4.xlsx"
file_DEGs_Batiuk = "D:/Documents/projects/SCZ_human_cortex/data/Gene_lists/Batiuk_et_al_2022/Supplementary_Dataset_Tables_1-4.xlsx" #path_main + "data/Batiuk_et_al_2022/Supplementary_Dataset_Tables_1-4.xlsx"
file_DEGs_Bast = path_project + "/5b_differential_gene_expression_analysis/output/DEGs/design_with_6_RUVs/T5_DEGs_per_cell_type.csv"#df_DESeq2_all_genes_all_CTs.csv"

alpha_vals = [0.05,0.1,0.3]
#min_subset_size = [1,5,10]

for i,alpha in enumerate(alpha_vals):
    #get p-values and log2FC of all datasets in one data frame
    DF_p, DF_l = ut.get_all_merged_DFs(opt_load_merged_DFs,path_results,file_DEGs_Ruzicka,file_DEGs_Batiuk,file_DEGs_Bast,alpha,opt_subset_sign_BAST_genes)
    
    
    #plot upset plot, merge CTs in classes for all data sets, use specific alpha cutoff
    #ut.plot_upset_gene_overlap_per_class(DF_p,alpha,path_results,min_subset_size[i],True)

    #DF_p_log10 = DF_p.copy()
    #for col in DF_p.columns.tolist():
    #    if col!="gene":
    #        DF_p_log10[col] = np.log10(DF_p[col])

    #DF_p_spear,DF_l_spear = ut.get_corr_matrix(DF_p_log10,DF_l,"spearman")
    #DF_p_pear,DF_l_pear = ut.get_corr_matrix(DF_p_log10,DF_l,"pearson")
    DF_p_spear,DF_l_spear = ut.get_corr_matrix(DF_p,DF_l,"spearman")

    for sel in ["all","Exc","Inh","NN"]:
        column_sel = DF_p_spear.columns.tolist()
        rows_sel = DF_p_spear.index.tolist()
        if sel == "Exc":
            column_sel = [c for c in column_sel if c.startswith("Ruzicka_Ex-") or (c.startswith("Batiuk_") and "LAMP5" not in c and "SST" not in c and "VIP" not in c and "PVALB" not in c)] 
            rows_sel = [r for r in rows_sel if r.startswith("Bast_Exc")]
        elif sel == "Inh":
            column_sel = [c for c in column_sel if c.startswith("Ruzicka_In-") or "LAMP5" in c or "SST" in c or "VIP" in c or "PVALB" in c] 
            rows_sel = [r for r in rows_sel if r.startswith("Bast_Inh")]
        elif sel == "NN":
            column_sel = [c for c in column_sel if c.startswith("Ruzicka") and not c.startswith("Ruzicka_In-") and not c.startswith("Ruzicka_Ex")] #"Oli" in c or "Ast" in c or "Mic" in c or "OPC" in c or "Peri" in c or "Endo" in c] 
            rows_sel = [r for r in rows_sel if not r.startswith("Bast_Exc") and not r.startswith("Bast_Inh")]

        ut.plot_correlations_as_heatmap(DF_l_spear,"_spearman_corr_log2FC",path_results,opt_order,"log2FC",alpha,column_sel,rows_sel,sel)
        ut.plot_correlations_as_heatmap(DF_p_spear,"_spearman_corr_adj_p",path_results,opt_order,"adjusted p-value",alpha,column_sel,rows_sel,sel)
        #ut.plot_correlations_as_heatmap(DF_l_pear,"_pearson_corr_log2FC",path_results,opt_order,"log2FC",alpha)


    
