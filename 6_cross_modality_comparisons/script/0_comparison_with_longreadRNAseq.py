# author: Lisa Bast
# date: 2023-09-11,  09:23:28
# version: 0.0.1
# about: compare DEGs with differential isoform genes from longread RNAseq (Abrantes et al)
#        Fig. S13

import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#settings
n_cl=15
alpha_cutoffs = [0.05,0.1,0.2,0.3]

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()

#paths:
path_project_main = path_code.replace("script","")
path_main = path_code.replace("6_cross_and_inter_modality_comparisons\\script","")
path_DEG_result = path_main + "/5b_differential_gene_expression_analysis/output/DEGs/16_CTs/design_with_6_RUVs/"
DEG_filename = "df_DESeq2_all_genes_all_CTs.csv"
path_genelists = path_main + "/5b_differential_gene_expression_analysis/data/"
gene_matrix_file = path_genelists + "geneMatrix_v2.tsv"
path_longread = path_genelists + "longreadRNAseq/"
pathways_DIUGs_filename = path_longread + "/go_hypergeoPathway_filtE_fdr10.tsv"
DIUG_filename = "sr_mm_PC12_filtE_gene_FDR.tsv"
path_filtered_data = path_main + "/4_data_integration_and_cell_type_annotation/output/"
path_filtered_data_per_CT = path_filtered_data + "/CT_specific_files/16_CTs/"

#load helper functions:
sys.path.append(path_code)
import utils as ut1
sys.path.append(path_main + "/5b_differential_gene_expression_analysis/code/")
import utils as ut2

#plot stats
ut1.plot_DIUG_DEG_stats(DEG_filename, DIUG_filename,path_longread,alpha_cutoffs,"without_tested")
    
#keep only sign
#DF = DF[np.logical_or(DF["bool_DIG"]==True,DF["padj"])]
#list common genes

#pie charts with CT info
for alpha_cutoff in alpha_cutoffs:
    #investigate results
    tested_genes_DIG_and_DEG,tested_genes_DEG_only,tested_genes_DIG_only,sign_genes_DIG_and_DEG,sign_genes_DEG_only,sign_genes_DIG_only = ut1.get_gene_DFs_for_alpha_cutoffs_longread(DEG_filename, DIUG_filename, alpha_cutoff, path_longread)
    #read DEG result:   
    DF_DEGs = pd.read_csv(path_DEG_result+'df_DESeq2_all_genes_all_CTs.csv')
    DF_DEGs = ut1.get_short_CT_names(DF_DEGs,"celltype")
    DF_DEGs["Color_sign"][np.logical_and(DF_DEGs['log2FoldChange']<0,DF_DEGs['padj']<=alpha_cutoff)]="darkred"
    DF_DEGs["Color_sign"][np.logical_and(DF_DEGs['log2FoldChange']>0,DF_DEGs['padj']<=alpha_cutoff)]="teal"
    
    if alpha_cutoff>0.05:
        DF_DEGs["longreadRNAseq"] = DF_DEGs["Gene_short"].isin(sign_genes_DIG_and_DEG)
        for ct in DF_DEGs["celltype"].unique():
            DF_DEGs_CT = DF_DEGs[DF_DEGs["celltype"]==ct]
            DF_DEGs_CT.set_index("Gene_short",inplace=True)
            #volcano plot with DIUGs as annotations
            ut1.plot_volcano(DF_DEGs_CT,ct,path_longread+"DEG_comparison/",alpha_cutoff,True,annotations="longreadRNAseq")
        #plot for sign DIUGs some stats per cell type
        DF_sel = DF_DEGs[DF_DEGs["longreadRNAseq"]==True]

        DF_num = ut2.get_number_of_deregulated_genes(DF_sel,"log2FoldChange","celltype_short","Gene_short",n_cl,path_filtered_data)
        DF_num_detailed = ut2.get_gene_type_for_deregulated_genes(DF_sel[["log2FoldChange","celltype_short","Gene_short"]], "log2FoldChange", gene_matrix_file, "celltype_short","Gene_short")
        DF_num_detailed = DF_num_detailed.sort_values(by="total",ascending=False)
        ut2.plot_detailed_stacked_bar(DF_num_detailed,True,path_longread+"DEG_comparison/",'_num_dereg_genes_'+str(alpha_cutoff).replace('.','_')+'.pdf')

        #plot heatmap
        #genes vs CTs log2FC
        ut1.plot_heatmap_log2FC_longread_comparison(DF_DEGs,path_longread+"DEG_comparison/figures/",alpha_cutoff,n_cl,path_filtered_data)

    if len(sign_genes_DIG_and_DEG):
        DF = DF_DEGs[DF_DEGs["Gene_short"].isin(sign_genes_DIG_and_DEG)].copy()
        DF['regulation_DEG'] = "not de-regulated"
        DF['regulation_DEG'][np.logical_and(DF['log2FoldChange']>0,DF['padj']<=alpha_cutoff)] = "up"
        DF['regulation_DEG'][np.logical_and(DF['log2FoldChange']<0,DF['padj']<=alpha_cutoff)] = "down"
        DF.drop(DF[DF['regulation_DEG']=="not de-regulated"].index, inplace=True)
        #fp.plot_number_sign_genes_overlapping_upset(DF,True,path_longread,alpha_cutoff,'longreadRNAseq',1)
        ut1.plot_number_sign_genes_overlapping_bar(DF,True,path_longread,'longreadRNAseq',alpha_cutoff)


#compare longread RNA-seq pathways with DEG pathways:
df_sn = pd.read_csv(path_DEG_result+DEG_filename)
CT_columns = df_sn.columns.tolist()
CT_columns = [c for c in CT_columns if c not in ["group","subgroup","geneset"]]

df_lr = pd.read_table(pathways_DIUGs_filename,index_col=0)
df_lr_sign = df_lr[df_lr["BH"]<0.05]

df_overlap = pd.merge(left = df_lr_sign[["geneset","BH"]], right = df_sn, on="geneset",how="left")
df_overlap["Number of sign celltypes"] = (df_overlap[CT_columns]<0.05).sum(axis=1)
df_overlap = df_overlap.sort_values(by = "Number of sign celltypes",ascending=False)
df_overlap.rename(columns={"BH":"longread_RNAseq_BH"},inplace=True)

# To DO: make sure DF is in correct format!
#heatmap of overlapping longread and snRNAseq pathways
df_overlap_with_MI = df_overlap[df_overlap.columns[df_overlap.columns!="longread_RNAseq_BH"]]
#sort by subgroup
df_overlap_with_MI.sort_values(by='group',inplace=True)
df_overlap_with_MI.set_index(['subgroup','geneset'],inplace=True)
#To do: fix
#df_overlap_with_MI = df_overlap_with_MI.fillna(0)
#get short ct names first:
df_overlap_with_MI = ut2.get_short_CT_names_as_columns(df_overlap_with_MI)
#sort columns
#CT_ordered = fg.get_CT_order(n_cl,path_filtered_data,True)
df_overlap_with_MI.drop(columns=["group"],inplace=True)

ut2.plot_triangular_heatmap_GSEA_pathways(df_overlap_with_MI,"GOs",True, path_longread+'pathway_comparison/', False, False, True, 15, path_filtered_data, title_str = 'significant longread RNAseq pathways', background_str = '')

#plot to how many significant longread pathways each cell type is mapping to
DF_CT = (df_overlap[CT_columns]<=0.05).sum().to_frame(name = "Number of significant longread RNAseq pathways")
DF_CT["CT_long"] = DF_CT.index.tolist()
DF_CT[["Celltype","Regulation"]] = DF_CT.CT_long.str.rsplit(pat='_',n=1,expand=True)
DF_CT = DF_CT.reset_index()
DF_CT.drop(columns=["CT_long","index"],inplace=True)
DF_CT["Celltype"] = DF_CT["Celltype"].str.replace("_"," ")

ut1.plot_common_pathway_stats_per_CT(DF_CT, "longread RNAseq",path_longread+'pathway_comparison/')
