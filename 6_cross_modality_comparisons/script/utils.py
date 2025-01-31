# author: Lisa Bast
# date: 2024-10-29,  14:17:55
# version: 0.0.1
# about: helper functions for 6_cross_and_inter_modality_comparisons

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from sctriangulate.colors import build_custom_continuous_cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

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

def get_gene_DFs_for_alpha_cutoffs_longread(DEG_filename, DIUG_filename, alpha_cutoff, path_longread):
    #DEGs
    df_DEGs = pd.read_csv(DEG_filename)
    DEG_tested_genes = np.unique(df_DEGs["Gene_short"].to_list())

    #DIGs
    df_DIG = get_DIUG_dataframe(DIUG_filename,alpha_cutoff)

    #number of genes tesed DIG only/ DEG only/ both
    DIG_tested_genes = df_DIG["Gene_short"].to_list()

    tested_genes_DIG_and_DEG = np.intersect1d(DIG_tested_genes,DEG_tested_genes)
    tested_genes_DIG_only = list(set(DIG_tested_genes) - set(tested_genes_DIG_and_DEG))
    tested_genes_DEG_only = list(set(DEG_tested_genes) - set(tested_genes_DIG_and_DEG))
    #merge data frames
    DF = pd.merge(left=df_DEGs, right= df_DIG, on = "Gene_short")
    #list genes sign
    DIG_genes = DF[DF["bool_DIG"]==True]["Gene_short"].tolist()
    DEG_genes = np.unique(DF[DF["padj"]<=alpha_cutoff]["Gene_short"].tolist())

    sign_genes_DIG_and_DEG = np.intersect1d(DIG_genes,DEG_genes)
    sign_genes_DEG_only = list(set(DEG_genes) - set(sign_genes_DIG_and_DEG))
    sign_genes_DIG_only = list(set(DIG_genes) - set(sign_genes_DIG_and_DEG))

    save_stats_DIUG_DEG(tested_genes_DIG_and_DEG,"tested_genes_DIUG_and_DEG_"+str(alpha_cutoff).replace('.','_'),path_longread)
    save_stats_DIUG_DEG(tested_genes_DIG_only,"tested_genes_DIUG_only_"+str(alpha_cutoff).replace('.','_'),path_longread)
    save_stats_DIUG_DEG(tested_genes_DEG_only,"tested_genes_DEG_only_"+str(alpha_cutoff).replace('.','_'),path_longread)
    save_stats_DIUG_DEG(sign_genes_DIG_and_DEG,"sign_genes_DIUG_and_DEG_"+str(alpha_cutoff).replace('.','_'),path_longread)
    save_stats_DIUG_DEG(sign_genes_DEG_only,"sign_genes_DEG_only_"+str(alpha_cutoff).replace('.','_'),path_longread)
    save_stats_DIUG_DEG(sign_genes_DIG_only,"sign_genes_DIUG_only_"+str(alpha_cutoff).replace('.','_'),path_longread)

    return tested_genes_DIG_and_DEG,tested_genes_DEG_only,tested_genes_DIG_only,sign_genes_DIG_and_DEG,sign_genes_DEG_only,sign_genes_DIG_only

def change_columnnames_in_df(df,add_str):
    for c in df.columns:
        df = df.rename(columns={c:add_str+c})
    return df

def shorten_gene_names(df):
    df.reset_index(inplace=True)
    df[['gene','ensgid']]  = df["Gene"].str.split("_",expand=True)
    df.drop(columns=["Gene","ensgid"],inplace=True)
    return df

def dec_to_rgb(rgb_tuple_dec):
    r = rgb_tuple_dec[0]*255.0
    b = rgb_tuple_dec[1]*255.0
    g = rgb_tuple_dec[2]*255.0
    return (r,b,g)
    #https://stackoverflow.com/questions/25404998/how-can-i-convert-an-rgb-input-into-a-decimal-color-code
    #return (r << 16) + (g << 8) + b

def get_DIUG_dataframe(DIUG_filename,alpha_cutoff):
    df_DIG = pd.read_table(DIUG_filename,index_col=0)
    df_DIG["bool_DIG_tested"] = True
    df_DIG["bool_DIG"] = False
    df_DIG.loc[df_DIG["FDR"]<=alpha_cutoff,"bool_DIG"] = True
    df_DIG.rename(columns={"associated_gene": "Gene_short"},inplace=True)

    return df_DIG

#save stats:
def save_stats_DIUG_DEG(list_genes,filename_str,path_longread):
    # add '\n' after each item of a list
    n_names = ["{}\n".format(i) for i in list_genes]
    with open(path_longread+'DEG_comparison/'+filename_str+'.txt', 'w') as f:
        f.writelines(n_names)

def get_CT_colors_as_df(n_clusters,path_filtered_data):
    if path_filtered_data.endswith('cellranger/'):
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
    else:
        print('This is not implemented yet.')
        DF = []
    return DF

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

def get_DEP_DEG_module_dataframe(DF,path_genelists,path_module_genes,path_DEG_result, DEG_result_filename,path_results_proteomics,filename_proteomics,alpha_val_DEPs,path_results_proteomics_hgnc,hgnc_file,opt_with_MT_genes,bordeaux_modules):
    DF_p = get_DEP_results(path_results_proteomics,filename_proteomics,1 ,"","", alpha_val_DEPs ,path_results_proteomics_hgnc+hgnc_file)

    #for the bordeaux modules: plot log2FC in proteomics and highlight MT genes and nuclear genes
    if opt_with_MT_genes==True:
        DF_BM = DF[DF["module_name"].isin(bordeaux_modules)]
    else:
        DF_BM =DF.copy()

    DF_BM_p = pd.merge(left=DF_BM, right=DF_p, left_on="accession",right_on="Ensembl gene ID",how="left")
    DF_BM_p.dropna(inplace=True) # not all genes in bordeaux modules map to proteins in DEP data set --> 89 unique genes
    DF_BM_p["MT_gene_status"]=DF_BM_p["Gene"].str.startswith("MT-")
    #DF_BM_p = fg.get_short_CT_names(DF_BM_p,"module_name")
    #DF_BM_p["module_name_short"] = DF_BM_p["module_name_short"].str.replace("_module_", " ").str.replace("_"," ")
    DF_BM_p = DF_BM_p.sort_values("log2FoldChange_DEP")
    DF_BM_p_plot = DF_BM_p[["gene_ids","log2FoldChange_DEP","MT_gene_status"]].drop_duplicates()
    DF_BM_p_plot.reset_index(inplace=True,drop=True)
    DF_BM_p_plot ['x1'] = DF_BM_p_plot .index
    
    #log2FC transcriptomics vs log2FC proteomics#
    DF_t = pd.read_csv(path_DEG_result+DEG_result_filename)
    DF_t.rename(columns={"log2FoldChange":"log2FoldChange_DEG","padj":"adjusted p-value_DEG"},inplace=True)
    DF_t["signed adjusted p-value_DEG"] = DF_t["adjusted p-value_DEG"]
    DF_t.loc[DF_t["log2FoldChange_DEG"]<0,"signed adjusted p-value_DEG"] = DF_t["signed adjusted p-value_DEG"]*(-1)
    DF_BM_t = pd.merge(left=DF_BM, right=DF_t, left_on=["celltype","gene_ids"],right_on=["celltype","Gene_short"])
    #ranked log2FC 
    DF_BM_t = DF_BM_t.sort_values("log2FoldChange_DEG")
    DF_BM_t.reset_index(inplace=True,drop=True)
    DF_BM_t['x1'] = DF_BM_t.index
    DF_BM_t["MT_gene_status"]=DF_BM_t["Gene_short"].str.startswith("MT-")
    
    DF_BM_p.rename(columns={"qvalue_DEP":"adjusted p-value_DEP"},inplace=True)
    DF_BM_t_p = pd.merge(left=DF_BM_t,right=DF_BM_p,on=["gene_ids","module_name","MT_gene_status"])
    DF_BM_t_p["signed adjusted p-value_DEP"] = DF_BM_t_p["adjusted p-value_DEP"]
    DF_BM_t_p.loc[DF_BM_t_p["log2FoldChange_DEP"]<0,"signed adjusted p-value_DEP"] = DF_BM_t_p["signed adjusted p-value_DEP"]*(-1)

    DF_p = DF_p[DF_p["qvalue_DEP"]<alpha_val_DEPs]
    DF_p["Gene"].to_csv(path_genelists+"sign_DEP_genes.csv",header=False, index=False)
    DF_p["Ensembl gene ID"].to_csv(path_genelists+"sign_DEP_accessions.csv",header=False, index=False)
    DF = pd.merge(left = DF_p[["Gene","log2FoldChange_DEP"]], right = DF,right_on="gene_ids",left_on="Gene",how="left")
    DF.dropna(subset=["module_name"],inplace=True)
    
    return DF_BM_p_plot,DF_BM_t,DF_BM_t_p,DF

def get_matrix_of_modules_vs_pathways(DF,opt_filter_genes,mode):
    DF["count"] = 1
    DF = get_short_CT_names(DF,"module_name")
    DF["module_name_short"] = DF["module_name_short"].str.replace("_module_"," ")
    DF["module_name_short"] = DF["module_name_short"].str.replace("_"," ")
    if opt_filter_genes=="no_filtering":
        M_df = DF.pivot(columns=["module_name_short"],index='accession', values=mode)
    else:
        M_df= DF.pivot(columns=["module_name_short"],index='gene_ids', values=mode)
    M_df[M_df.isna()]=0
    if mode=="count":
        opt_binary_plot=True
    else:
        opt_binary_plot = False
    return M_df,opt_binary_plot

def get_MT_genes_from_DEG_and_DAP(path_DEG_result,DEG_filename,path_results_proteomics_DAPs,filename_proteomics_DAPs,hgnc_file,path_filtered_data):

    #read_DEGs, extract mitochondrial DNA gene rows:
    DEGs = pd.read_csv(path_DEG_result + DEG_filename)
    #for MT genes: rank by log2FC and plot sorted log2FC
    DEG_MT = DEGs[DEGs["Gene_short"].str.startswith("MT-")]
    #add CT colors:
    #DF_colors = fg.get_CT_colors_as_df(n_cluster,path_colors)
    DEG_MT["celltype"] = DEG_MT["celltype"].str.replace("_"," ")
    DEG_MT = get_short_CT_names(DEG_MT,"celltype")
    #add proteomics
    DF_p_sel = get_DEP_results(path_results_proteomics_DAPs,filename_proteomics_DAPs,1,[''],[''],1,hgnc_file)
    DEP_MT = DF_p_sel[DF_p_sel["Gene"].str.startswith("MT-")]
    DEP_MT.rename(columns={"Gene":"Gene_short","log2FoldChange_DEP":"log2FoldChange","qvalue_DEP":"padj"}, inplace=True)
    DEP_MT["Color_sign"] = "grey"
    DEP_MT.loc[np.logical_and(DEP_MT["padj"]<=0.05,DEP_MT["log2FoldChange"]<0),"Color_sign"]="teal"
    DEP_MT.loc[np.logical_and(DEP_MT["padj"]<=0.05,DEP_MT["log2FoldChange"]>0),"Color_sign"]="darkred"
    DEP_MT["celltype_short"]="Proteomics"
    DEG_MT = pd.merge(DEG_MT,DEP_MT,on=["Gene_short","log2FoldChange","celltype_short","padj","Color_sign"],how="outer")
    #change order of CTs
    CT_ordered = get_CT_order(15,path_filtered_data,True)
    CT_ordered = ["Proteomics"]+CT_ordered
    #DEG_MT["celltype_short"] = pd.Categorical(DEG_MT.celltype_short, ordered=True, categories=CT_ordered)
    DEG_MT = DEG_MT.set_index("celltype_short")
    DEG_MT = DEG_MT.reindex(index = CT_ordered)
    DEG_MT = DEG_MT.reset_index()

    return DEG_MT


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

def plot_number_sign_genes_overlapping_bar(DF,opt_save,path_results,input_data,alpha_val):
    DF['up-regulated']=False
    DF['up-regulated'][DF["regulation_DEG"]=="up"]=True
    DF['down-regulated']=False
    DF['down-regulated'][DF["regulation_DEG"]=="down"]=True
    DF['total']=True
    DF_number = DF.groupby('celltype_short').sum()[['up-regulated','down-regulated','total']].sort_values(by='total',ascending=False)
    DF_number.index = DF_number.index.str.replace('_',' ')
    axes=DF_number[['down-regulated','up-regulated']].plot(kind='barh',stacked=True,color={"down-regulated":'darkred',"up-regulated":'teal'},edgecolor="white")

    #axes=DF_number[['down-regulated','up-regulated']].plot(kind='barh',stacked=True,color={"down-regulated":'teal',"up-regulated":'darkred'},edgecolor="white")
    axes.ticklabel_format(style='plain', axis='x') 
    axes.set_ylabel("Number of genes")
    axes = remove_frame_keep_axes(axes)
    plt.gca().invert_yaxis()
    if opt_save:
        if input_data=="proteomics":
            path=path_results
        else:
            path=path_results+'DEG_comparison/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'bar_DEGs_CT_overlap_with_'+input_data+'_'+str(alpha_val).replace('.','_')+'.pdf', bbox_inches='tight')
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

def plot_DIUG_DEG_stats(DEG_filename, DIUG_filename,path_longread,alpha_cutoffs,opt):
    if opt!="with_tested":
        x=[]
        y1=[]
        y2=[]
        y3=[]
        n_cols=len(alpha_cutoffs)
    else:
        n_cols=len(alpha_cutoffs)+1
    #get data in right format:
    for i,alpha_cutoff in enumerate(alpha_cutoffs):
        tested_genes_DIG_and_DEG,tested_genes_DEG_only,tested_genes_DIG_only,sign_genes_DIG_and_DEG,sign_genes_DEG_only,sign_genes_DIG_only = get_gene_DFs_for_alpha_cutoffs_longread(DEG_filename, DIUG_filename, alpha_cutoff,path_longread)
        
        if i==0 and opt=="with_tested":
            y1 = np.array([len(tested_genes_DIG_and_DEG),len(sign_genes_DIG_and_DEG)])
            y2 = np.array([len(tested_genes_DEG_only),len(sign_genes_DEG_only)])
            y3 = np.array([len(tested_genes_DIG_only),len(sign_genes_DIG_only)])
            x= ["tested","significant (alpha="+str(alpha_cutoff)+")"]
        else:
            y1 = np.append(y1,len(sign_genes_DIG_and_DEG))
            y2 = np.append(y2,len(sign_genes_DEG_only))
            y3 = np.append(y3,len(sign_genes_DIG_only))
            x = x + ["significant (alpha="+str(alpha_cutoff)+")"]
    #plot:
    fig,ax = plt.subplots(nrows=1,ncols=n_cols,sharey=True)
    fig.tight_layout()
    for i in range(0,n_cols):
        ax[i].bar(x[i],y1[i],color="#003366",label="DIUG and DEG")
        ax[i].bar(x[i],y2[i],bottom = y1[i], color="#E0E0E0",label="DEG only")
        ax[i].bar(x[i],y3[i],bottom = y1[i]+y2[i], color="#A0A0A0",label="DIUG only")
        if i==0:
            ax[i].set_ylabel("Number of Genes")
        elif i==n_cols-1:
            ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax[i]=remove_frame_keep_axes(ax[i])
    fig.show()
    path=path_longread+'DEG_comparison/'
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    fig.savefig(path + 'DEG_DIUG_comparison_stats_'+opt+'.pdf',bbox_inches='tight')
    plt.close(fig)
    plt.clf() 

def plot_common_pathway_stats_per_CT(DF_CT, plot_for_str,results_path):
    DF_pivot = pd.pivot_table(DF_CT,index = "Celltype",columns="Regulation",values="Number of significant "+plot_for_str+" pathways")
    CT_names = DF_pivot.index.tolist()
    CT_names_reordered = np.array(CT_names)[np.char.startswith(CT_names,'Exc')].tolist() + np.array(CT_names)[np.char.startswith(CT_names,'Inh')].tolist() + np.array(CT_names)[np.logical_and(~np.char.startswith(CT_names,'Exc'),~np.char.startswith(CT_names,'Inh'))].tolist()
    DF_pivot = DF_pivot.reindex(CT_names_reordered)
    ax = DF_pivot.plot.bar(stacked=True,color={"up":"darkred","down":"teal"})
    ax.set_ylabel("Number of common pathways")
    plt.tight_layout()

    isExist = os.path.exists(results_path)
    if not isExist:
        os.makedirs(results_path)
    plt.savefig(results_path + 'pathways_'+plot_for_str+'_DEG_comparison_stats.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

def plot_heatmap_log2FC_longread_comparison(DF_DEGs,path_longread,alpha_cutoff,n_cl,path_filtered_data):
    DF_HM = DF_DEGs[DF_DEGs["longreadRNAseq"]==True][["Gene_short","log2FoldChange","celltype_short"]]
    DF_HM_p = pd.pivot_table(DF_HM, values="log2FoldChange", columns="celltype_short",index="Gene_short")
    CM = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('teal'))).tolist(),
                                      np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                      np.array(dec_to_rgb(colors.to_rgb('darkred'))).tolist())
    CM.set_bad("lightgrey")
    CTs_ordered = get_CT_order(n_cl,path_filtered_data,True)
    CTs_ordered_ = [ct.replace(' ','_') for ct in CTs_ordered]
    DF_HM_p = DF_HM_p[CTs_ordered_]
    plt = loadPltSettings(4,3)
    cl = sns.clustermap(data=DF_HM_p.fillna(value=0),col_cluster=False)
    genes_ordered = DF_HM_p.index[cl.dendrogram_row.reordered_ind].tolist()
    DF_HM_p_ordered = DF_HM_p.loc[genes_ordered]
    plt.close()
    max_abs_val = np.max(np.abs([DF_HM_p_ordered.max().max(),DF_HM_p_ordered.min().min()]))
    g=sns.heatmap(data=DF_HM_p_ordered.transpose(), center=0,cmap=CM,square=True,vmin=max_abs_val*(-1),vmax=max_abs_val)
    g.axes.set_xticklabels(g.axes.get_xmajorticklabels(), fontsize = 4)
    g.axes.set_yticklabels(g.axes.get_ymajorticklabels(), fontsize = 4)
    if not os.path.exists(path_longread):
        os.makedirs(path_longread)         
    plt.savefig(path_longread+'_heatmap_log2F_dereg_genes_'+str(alpha_cutoff).replace('.','_')+'.pdf', bbox_inches='tight')
    plt.close()
    plt.clf() 

def plot_clustermap_module_pathway_comparison(df_overlap_p,opt_order_modules,path_results,compare_with_str,opt_test,opt_visualize_p_val):
    if opt_visualize_p_val:
        add_str2 = "_p_FDR_group"
        cmap_us = sns.color_palette("blend:#005675,#FFFFFF", as_cmap=True)
    else:
        add_str2 = ""
        cmap_us = sns.color_palette("blend:#FFFFFF,#005675", as_cmap=True)
    #square clustermap
    if compare_with_str=="proteomics":
        fw, fh = 10, 10
    else:
        fw, fh = 4, 4
    #remove modules with no pw overlap:
    df_overlap_p = df_overlap_p[df_overlap_p.columns[df_overlap_p.sum(axis=0)!=0]]
    if opt_order_modules:
        g = sns.clustermap(df_overlap_p, cmap=cmap_us,figsize=(6, 6),yticklabels=1,xticklabels=1,col_cluster = False)
        add_str = "_ordered"
    else:
        g = sns.clustermap(df_overlap_p, cmap=cmap_us,figsize=(6, 6),yticklabels=1,xticklabels=1)
        add_str = ""
    g.cax.set_visible(False)
    dx, dy = np.diff(g.ax_heatmap.transData.transform([[0, 1], [1, 0]]), axis=0).squeeze()
    if dx > dy:
        fh = fh * dx / dy
    else:
        fw = fw * dy / dx
    g.fig.set_size_inches(fw, fh)
    dx, dy = np.diff(g.ax_heatmap.transData.transform([[0, 1], [1, 0]]), axis=0).squeeze()
    g.fig.savefig(path_results+"Clustermap_"+compare_with_str+"_pathways_per_module_all_GOs"+'_'+opt_test+add_str+add_str2+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

def plot_clustermaps(DF_C,DF_P,M_df,path_module_genes,opt_test,opt_filter_genes,n_cl,path_filtered_loom_data,gene_matrix_file,opt_binary_plot):
    if opt_binary_plot:
        color_map = "Blues"
        add_str = "_binary"
    else:
        color_map = sns.palplot(sns.diverging_palette(220, 20, n=7))
        add_str = "_log2FC"

    sns.clustermap(data=DF_C, annot=None, cmap=color_map,yticklabels=True,xticklabels=True)
    plt.savefig(path_module_genes+'Correlation_of_Module_wrt_genes_they_contain_'+opt_test+'_'+opt_filter_genes+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()
    
    cg = sns.clustermap(data=DF_P, annot=None, cmap='Blues',vmin = 0, vmax=100, yticklabels=True,xticklabels=True)#, cbar_kws={"ticks":[0,100]})
    #enlarge figure
    cg.fig.set_size_inches(8,10)
    # make some space to the right in the figure
    cg.gs.update(top=0.95)
    # divide existing axes
    divider = make_axes_locatable(cg.ax_heatmap)
    divider2 = make_axes_locatable(cg.ax_row_dendrogram)
    # create new axes for bar plot 
    ax = divider.append_axes("top", size="30%", pad=0.1)
    # create empty space of same size as bar plot axes (don't use this space)
    nax = divider2.new_vertical(size="20%", pad=1.7)

    # Sort the values for the bar plot to have the same order as clusters
    pd.options.display.float_format = '{:.2f}'.format
    colors = get_CT_colors_as_df(n_cl,path_filtered_loom_data)
    colors = get_short_CT_names(colors,"celltype")
    colors["celltype_short"] = colors["celltype_short"].str.replace("_"," ")

    M_df.reset_index(inplace=True)
    #Number_genes_with_type = fg.get_gene_type_for_deregulated_genes(Number_genes, [], gene_matrix_file,[])
    M_df[["celltype_short","module"]] = M_df["module_name_short"].str.rsplit(pat=" ",n=1,expand=True)
    M_df = pd.merge(left = M_df, right=colors,on="celltype_short")
    #sort cts:
    #M_df["module_name_short"] =
    M_df.set_index(keys="module_name_short",drop=True, inplace=True)
    M_df.drop(columns=["celltype","celltype_short","module"],inplace=True)
    #load order of modules:
    module_list_ordered = DF_P.index[cg.dendrogram_row.reordered_ind].tolist()

    M_df_ordered = M_df.loc[module_list_ordered]
    if opt_filter_genes=="no_filtering":
        genes_columns = M_df_ordered.columns.tolist()
        genes_columns.remove("color")
        #determine number of lncRNA and protein coding genes
        #load geneMatrix file
        GM = pd.read_csv(gene_matrix_file,sep='\t')
        lncRNA_genes_GM = GM[GM["gene_type"]=="lncRNA"]["ensgid"].tolist()
        pc_genes_GM = GM[GM["gene_type"]=="protein_coding"]["ensgid"].tolist()
        lnc_genes = [g for g in genes_columns if g in lncRNA_genes_GM]
        pc_genes = [g for g in genes_columns if g in pc_genes_GM]
        other_genes = [g for g in genes_columns if g not in lnc_genes and g not in pc_genes]
        #GM[GM["gene_name"].isin(other_genes)]["gene_type"]
        #merge genetype info from GM into DF_DEG_result_red
        Number_genes = M_df_ordered.loc[:,lnc_genes].sum(axis=1).to_frame(name="lncRNA").merge(M_df_ordered.loc[:,pc_genes].sum(axis=1).to_frame(name="protein coding"),how="left", left_on="module_name_short",right_on="module_name_short")
        Number_genes = Number_genes.merge(M_df_ordered.loc[:,other_genes].sum(axis=1).to_frame(name="other"),how="left", left_on="module_name_short",right_on="module_name_short")
        Number_genes = Number_genes.merge(M_df_ordered.loc[:,"color"].to_frame(name="color"),how="left", left_on="module_name_short",right_on="module_name_short")
        Number_genes["color"].tolist()
        Number_genes[["lncRNA","protein coding", "other"]].plot(kind='bar',
                                                                stacked=True,
                                                                xlabel='Number of genes',
                                                                ax=ax,
                                                                color = {"lncRNA":"lightgrey","protein coding":"darkgrey", "other":"black"},
                                                                #color={'down-regulated protein coding':'teal','down-regulated lncRNA':'#104E4C','up-regulated protein coding':'darkred','up-regulated lncRNA':'#5B0B0B'},
                                                                edgecolor="white") 
    else:
        Number_genes = M_df_ordered.loc[:,M_df_ordered.columns != "color"].sum(axis=1)
        #pl = M_df_ordered.loc[:,M_df_ordered.columns != "color"].sum(axis=1).plot(kind="barh",color=M_df_ordered["color"],xlabel = "Number of Genes")
        #pl.ticklabel_format(axis="x", style="plain")
        #pl = remove_frame_keep_axes(pl)
        
        #target = [t.get_text() for t in np.array(cg.ax_heatmap.get_xticklabels())]
        #ind= np.array([list(M_df.index.values).index(t) for t in target])

        # plot bar plot in ax
        ax.bar(np.arange(len(module_list_ordered)), Number_genes,color = M_df_ordered["color"])
    ax.set_xticklabels([])
    ax = remove_frame_keep_axes(ax)
    ax.set_xlim(-0.5,len(Number_genes.index)-.5)
    #ax.invert_yaxis()

    plt.savefig(path_module_genes+'Module_genes_overlap_in_percent_'+opt_test+'_'+opt_filter_genes+add_str+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

    return module_list_ordered

def plot_integrated_analysis_maps(DF_pw_pivot,GO_str,results_path):
    CM = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('lightgrey'))).tolist(),
                                      np.array(dec_to_rgb(colors.to_rgb('darkblue'))).tolist())
    CM.set_bad("white")
    sns.heatmap(DF_pw_pivot,cmap=CM, cbar_kws={'label': "Proportion pathways per analysis"})
    plt.savefig(results_path+'heatmap_number_pathways_'+GO_str+'.pdf')
    CM = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('lightgrey'))).tolist(),
                                         np.array(dec_to_rgb(colors.to_rgb("white"))).tolist(),
                                         np.array(dec_to_rgb(colors.to_rgb('darkblue'))).tolist())
    DF_ratio = DF_pw_pivot/DF_pw_pivot.sum(axis=0)
    sns.clustermap(DF_ratio.fillna(0),cmap=CM, center= DF_ratio.min().min(), cbar_kws={'label': "Proportion pathways per analysis"},row_cluster=True, col_cluster=False, square=True, yticklabels=True)
    plt.savefig(results_path +'clustermap_ratio_pathways'+GO_str+'.pdf')
    plt.close()
    plt.clf() 

def plot_clustermap_modules_vs_genes(M_df,path_module_genes,opt_test,opt_filter_genes,opt_binary_plot):
    if opt_filter_genes == "no_filtering":
        bool_cluster_rows = False
    else:
        bool_cluster_rows = True
    
    if opt_binary_plot:
        color_map = "Blues"
        add_str = "_binary"
    else:
        color_map = sns.diverging_palette(220, 20, n=100,as_cmap=True)
        add_str = "_log2FC"
    #to do: fp.plot_number_genes_per_module(M_df,n_cl,path_module_genes,path_data,opt_test,opt_filter_genes) --> add to x-axis instead of clustermap
    #for opt_filter_genes in ["DEP_genes","DIUGs"]: take order of modules instead of clustering --> plot heatmap
    cluster = sns.clustermap(data=M_df, 
                            annot=None, 
                            cmap=color_map,
                            center=0.00,
                            yticklabels=True, 
                            xticklabels=True, 
                            col_cluster = True,
                            row_cluster = bool_cluster_rows,
                            figsize=(10,8))

    plt.savefig(path_module_genes+'Clustermap_modules_vs_genes_'+opt_test+'_'+opt_filter_genes+add_str+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf() 

def plot_MT_genes_log2FC_per_CT(DEG_MT,path_results):
    fig,ax = plt.subplots(2,1,figsize=(15,5),sharex=False)
    pal = {"grey":"grey","teal" :"teal", "darkred":"darkred"}
    sns.stripplot(data=DEG_MT,x="Gene_short",y="log2FoldChange",hue = "Color_sign",palette=pal,ax=ax[0])
    ax[0].axhline(y=0,linestyle="--",color="grey") 
    sns.stripplot(data=DEG_MT,x="celltype_short",y="log2FoldChange",hue = "Color_sign",palette=pal,ax=ax[1])
    ax[1].tick_params(axis="x",rotation=90)
    ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax[1].axhline(y=0,linestyle="--",color="grey") 
    # plot the mean line
    sns.boxplot(showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 2},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10,
                x="celltype_short",
                y="log2FoldChange",
                data=DEG_MT,
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=ax[1])
    fig.savefig(path_results+"MT_genes_log2FC_per_CT.pdf",bbox_inches='tight')
    plt.close()
    plt.clf()  
    
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