# author: Lisa Bast
# date: 2025-03-31,  14:49:21
# version: 0.0.1
# about: helper functions for comparison of DEGs in Ruzika and Batiuk papers with our DEGs

import os
import pandas as pd
import seaborn as sns
from sctriangulate.colors import build_custom_continuous_cmap
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import numpy as np
import math
from upsetplot import plot, from_indicators

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

def plot_correlations_as_heatmap(df,add_str,results_path,opt_order,opt_value,alpha,column_sel,row_sel,sel):

    #drop columns and rows that are not in col_sel
    df = df.loc[row_sel,column_sel].copy()

    #clustermaps of correlation matrixes:
    CM = build_custom_continuous_cmap(np.array(dec_to_rgb(colors.to_rgb('palevioletred'))).tolist(),
                                        np.array(dec_to_rgb(colors.to_rgb('white'))).tolist(),
                                        np.array(dec_to_rgb(colors.to_rgb('mediumaquamarine'))).tolist())
    CM.set_bad("lightgrey")
    if opt_order == "cluster_map":
        g = sns.clustermap(df.fillna(0),
                            vmin=-1,
                            vmax=1,
                            center=0.00,
                            yticklabels=True, 
                            xticklabels=True, 
                            col_cluster = True,
                            row_cluster = True)
        CT_list_ordered = df.index[g.dendrogram_row.reordered_ind].tolist()
        df = df.reindex(CT_list_ordered)
        df = df[CT_list_ordered]
        plt.close()
    sns.set(font_scale=0.3)
    if sel=="NN":
        fig, ax = plt.subplots(figsize=(2,10)) 
        ann_size=4
    elif sel=="all":
        fig, ax = plt.subplots(figsize=(10,10)) 
        ann_size=2
    else:
        fig, ax = plt.subplots(figsize=(5,10)) 
        ann_size=3
    ma = df.max().max()
    mi = df.min().min()
    ax = sns.heatmap(df,
                vmin = (-1)*np.max(np.abs([mi,ma])),
                vmax = np.max(np.abs([mi,ma])),
                center = 0.00,
                yticklabels=True, 
                xticklabels=True, 
                square=True,
                cmap=CM,
                cbar_kws={"label":"correlation of " + opt_value + "s"},
                annot=True,
                annot_kws={"fontsize":ann_size},
                ax=ax)
    
    #add correlation value to cell if in 5th or 95th percentile
    for t in ax.texts:
        TH_95_percentile = np.nanpercentile(np.concatenate(df.to_numpy()),95)
        if float(t.get_text())>=TH_95_percentile or float(t.get_text())<(-1)*TH_95_percentile:
            t.set_text(str(np.round(float(t.get_text()),2))) #if the value is greater than 95th percentile or below -95th percentile then I set the text 
        else:
            t.set_text("") # if not it sets an empty text

    for m, tick_label in enumerate(ax.axes.get_yticklabels()):
        if tick_label.get_text().startswith("Bast") and tick_label.get_text().endswith(" (top 30 genes)"):
            tick_label.set_color("grey")
        else:
            tick_label.set_color("black")
        y_max = m

    #add vertical line between two data sets:
    j1 = np.nan
    j2 = np.nan
    j3 = np.nan
    for i, tick_label in enumerate(ax.axes.get_xticklabels()):
        if tick_label.get_text().startswith("Bast"):
            tick_label.set_color("black")
        elif tick_label.get_text().startswith("Batiuk_"):
            tick_label.set_color("black")
            j1=i
        elif tick_label.get_text().startswith("Ruzicka_"):
            tick_label.set_color("black")
            j2=i
        else:
            tick_label.set_color("black")
            j3=i
    lines_vec = [j1,j2,j3]
    lines_vec = [value for value in lines_vec if not math.isnan(value)]
    for n in range(0,len(lines_vec)):
        ax.axvline(x=np.sort(lines_vec)[n]+1, ymin = 0, ymax = y_max, linewidth=1, color="black")
    isExist = os.path.exists(results_path)
    if not isExist:
        os.makedirs(results_path)
    plt.savefig(results_path+"heatmap_"+sel+"_"+opt_order+add_str+'_'+str(alpha).replace("0.","")+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

def read_dfs(file_DEGs_Ruzicka,file_DEGs_Batiuk,file_DEGs_Bast):
    #read data frames
    #(1) Batiuk
    #xls = pd.ExcelFile(file_DEGs_Batiuk)
    DF_Batiuk = pd.read_excel(file_DEGs_Batiuk,sheet_name = "Table 1",skiprows=2) 
    DF_Batiuk.rename(columns={"Gene":"gene"},inplace=True)
    #DF_Batiuk_res_2 = DF_Batiuk[["Cell_type.1","Gene.1","padj.1","log2FoldChange.1"]].copy()
    #DF_Batiuk_res_2.drop_duplicates(keep='first',inplace=True)
    #DF_Batiuk_res_2.dropna(axis=0,inplace=True)
    #DF_Batiuk_res_2.rename(columns={"Gene.1":"gene"},inplace=True)

    DF_Batiuk_p_res1 = pd.pivot(DF_Batiuk,columns="Cell_type",index="gene",values="padj")
    #DF_Batiuk_p_res2 = pd.pivot(DF_Batiuk_res_2,columns="Cell_type.1",index="gene",values="padj.1")
    DF_Batiuk_p_res1 = change_columnnames_in_df(DF_Batiuk_p_res1,"Batiuk_")
    #DF_Batiuk_p_res2 = change_columnnames_in_df(DF_Batiuk_p_res2,"Batiuk_res2_")

    DF_Batiuk_l_res1 = pd.pivot(DF_Batiuk,columns="Cell_type",index="gene",values="log2FoldChange")
    #DF_Batiuk_l_res2 = pd.pivot(DF_Batiuk_res_2,columns="Cell_type.1",index="gene",values="log2FoldChange.1")
    DF_Batiuk_l_res1 = change_columnnames_in_df(DF_Batiuk_l_res1,"Batiuk_")
    #DF_Batiuk_l_res2 = change_columnnames_in_df(DF_Batiuk_l_res2,"Batiuk_res2_")

    #(2) Ruzicka
    xls = pd.ExcelFile(file_DEGs_Ruzicka)
    CTs_Ruzicka = xls.sheet_names
    bool_first = True
    for ct in CTs_Ruzicka:
        DF_Ruzicka_tmp = pd.read_excel(file_DEGs_Ruzicka,sheet_name = ct)
        DF_Ruzicka_tmp["celltype_Ruzicka"] = ct
        if bool_first:
            DF_Ruzicka = DF_Ruzicka_tmp.copy()
            bool_first = False
        else:
            DF_Ruzicka = pd.concat([DF_Ruzicka,DF_Ruzicka_tmp])

    DF_Ruzicka_p_meta = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="Meta_adj.P.Val")
    DF_Ruzicka_p_MtSinai = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="MtSinai_p_adj.loc")
    DF_Ruzicka_p_McLean = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="McLean_p_adj.loc")
    DF_Ruzicka_p_meta = change_columnnames_in_df(DF_Ruzicka_p_meta,"Ruzicka_")
    DF_Ruzicka_p_MtSinai = change_columnnames_in_df(DF_Ruzicka_p_MtSinai,"Ruzicka_MtSinai_")
    DF_Ruzicka_p_McLean = change_columnnames_in_df(DF_Ruzicka_p_McLean,"Ruzicka_McLean_")
    DF_Ruzicka_p_meta .reset_index(inplace=True)
    DF_Ruzicka_p_MtSinai .reset_index(inplace=True)
    DF_Ruzicka_p_McLean.reset_index(inplace=True)

    DF_Ruzicka_l_meta = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="Meta_logFC")
    DF_Ruzicka_l_MtSinai = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="MtSinai_logFC")
    DF_Ruzicka_l_McLean = pd.pivot(DF_Ruzicka,columns="celltype_Ruzicka",index="gene", values="McLean_logFC")
    DF_Ruzicka_l_meta = change_columnnames_in_df(DF_Ruzicka_l_meta,"Ruzicka_")
    DF_Ruzicka_l_MtSinai = change_columnnames_in_df(DF_Ruzicka_l_MtSinai,"Ruzicka_MtSinai_")
    DF_Ruzicka_l_McLean = change_columnnames_in_df(DF_Ruzicka_l_McLean,"Ruzicka_McLean_")
    DF_Ruzicka_l_meta.reset_index(inplace=True)
    DF_Ruzicka_l_MtSinai.reset_index(inplace=True)
    DF_Ruzicka_l_McLean.reset_index(inplace=True)

    #(3) Bast:
    DF_Bast = pd.read_csv(file_DEGs_Bast)
    DF_Bast_p = pd.pivot(DF_Bast,columns="celltype", index="Gene", values="padj")
    DF_Bast_p = change_columnnames_in_df(DF_Bast_p,"Bast_")
    DF_Bast_p = shorten_gene_names(DF_Bast_p)

    DF_Bast_l = pd.pivot(DF_Bast,columns="celltype", index="Gene", values="log2FoldChange")
    DF_Bast_l = change_columnnames_in_df(DF_Bast_l,"Bast_")
    DF_Bast_l = shorten_gene_names(DF_Bast_l)

    return DF_Ruzicka_p_meta,DF_Ruzicka_p_MtSinai,DF_Ruzicka_p_McLean,DF_Ruzicka_l_meta,DF_Ruzicka_l_MtSinai,DF_Ruzicka_l_McLean,DF_Batiuk_p_res1,DF_Batiuk_l_res1,DF_Bast_p,DF_Bast_l

def get_all_merged_DFs(opt_load_merged_DFs,results_path,file_DEGs_Ruzicka,file_DEGs_Batiuk,file_DEGs_Bast, alpha,opt_subset_sign_BAST_genes):
    if opt_load_merged_DFs:
        DF_p_sel = pd.read_csv(results_path+"DF_p_"+str(alpha).replace("0.","")+".csv")
        DF_l_sel = pd.read_csv(results_path+"DF_l_"+str(alpha).replace("0.","")+".csv")
    else:
        DF_Ruzicka_p_meta,DF_Ruzicka_p_MtSinai,DF_Ruzicka_p_McLean,DF_Ruzicka_l_meta,DF_Ruzicka_l_MtSinai,DF_Ruzicka_l_McLean,DF_Batiuk_p_res1,DF_Batiuk_l_res1,DF_Bast_p,DF_Bast_l = read_dfs(file_DEGs_Ruzicka,file_DEGs_Batiuk,file_DEGs_Bast)
        
        ## merge all p and all l DFs together
        DF_p = pd.merge(left = DF_Bast_p, right = DF_Ruzicka_p_meta, how="outer",on="gene")
        DF_p = pd.merge(left = DF_p, right = DF_Ruzicka_p_MtSinai, how="outer",on="gene")
        DF_p = pd.merge(left = DF_p, right = DF_Ruzicka_p_McLean, how="outer",on="gene")
        DF_p = pd.merge(left = DF_p, right = DF_Batiuk_p_res1, how="outer",on="gene")
        
        DF_l = pd.merge(left = DF_Bast_l, right = DF_Ruzicka_l_meta, how="outer",on="gene")
        DF_l = pd.merge(left = DF_l, right = DF_Ruzicka_l_MtSinai, how="outer",on="gene")
        DF_l = pd.merge(left = DF_l, right = DF_Ruzicka_l_McLean, how="outer",on="gene")
        DF_l = pd.merge(left = DF_l, right = DF_Batiuk_l_res1, how="outer",on="gene")
        #DF_l = pd.merge(left = DF_l, right = DF_Batiuk_l_res2, how="outer",on="gene")

        #DF_p = pd.merge(left = DF_p, right = DF_Batiuk_p_res2, how="outer",on="gene")
        if opt_subset_sign_BAST_genes:
            Bast_cols = DF_l.columns[DF_l.columns.str.startswith("Bast_")]
            for col in Bast_cols:
                #genes that are not sign for this CT should be set to np.nan in both data frames
                genes_drop_CT = DF_p[DF_p[col]>alpha]["gene"].tolist()
                n_genes_kept = len(DF_p[DF_p[col]<=alpha]["gene"].tolist())
                if n_genes_kept<30:
                    #identify top 30 most significant genes
                    DF_p_sorted = DF_p[col].sort_values()
                    genes_include = DF_p.iloc[DF_p_sorted.index[0:30]]["gene"].tolist()
                    #remove those from genes_drop_CT
                    genes_drop_CT = [g for g in genes_drop_CT if g not in genes_include]
                #set value to nan if this gene and CT info is to be dropped:
                DF_p.loc[DF_p["gene"].isin(genes_drop_CT),col] = np.nan
                DF_l.loc[DF_l["gene"].isin(genes_drop_CT),col] = np.nan
                
                #indicate in column name how many genes were kept (to calculate correlation)
                if n_genes_kept<30:
                    DF_p.rename(columns={col:col+" (top 30 genes)"},inplace=True)
                    DF_l.rename(columns={col:col+" (top 30 genes)"},inplace=True)
                else:
                    DF_p.rename(columns={col:col+" ("+str(n_genes_kept)+" genes )"},inplace=True)
                    DF_l.rename(columns={col:col+" ("+str(n_genes_kept)+" genes )"},inplace=True)

        #DF_p_sign_bool = DF_p[[c for c in DF_p.columns.tolist() if c!="gene"]]<=alpha
        #genes_keep = DF_p["gene"][DF_p_sign_bool.sum(axis=1)>=1].tolist()

        #DF_p_sel = DF_p[DF_p["gene"].isin(genes_keep)]
        #DF_l_sel = DF_l[DF_l["gene"].isin(genes_keep)]
        #if all col except gene are nan --> drop row:

        #drop all rows (genes) for which all CTs are nan
        DF_p_sel = DF_p.dropna(subset=DF_p.columns[DF_p.columns.str.startswith("Bast_")],how="all")
        DF_l_sel = DF_l.dropna(subset=DF_l.columns[DF_l.columns.str.startswith("Bast_")],how="all")

        #export DF_p and DF_l
        isExist = os.path.exists(results_path)
        if not isExist:
            os.makedirs(results_path)  
        DF_p_sel.to_csv(results_path+"DF_p_"+str(alpha).replace("0.","")+".csv",index=False)
        DF_l_sel.to_csv(results_path+"DF_l_"+str(alpha).replace("0.","")+".csv",index=False)

    
    #make DF_p signed
    mask = DF_l_sel[[c for c in DF_p_sel.columns.tolist() if c!="gene"]]<0
    mask_sign = mask*(-1) + ~mask*(1)
    DF_p_sel_signed=DF_p_sel.copy()
    DF_p_sel_signed[[c for c in DF_p_sel.columns.tolist() if c!="gene"]] = DF_p_sel_signed[[c for c in DF_p_sel.columns.tolist() if c!="gene"]]*mask_sign

    return DF_p_sel_signed,DF_l_sel

def make_DF_concise(DF):
    import re
    #drop rows that are not Bast et al:
    id_list = DF.index.tolist()
    row_ids_keep = [rowname for rowname in id_list if (rowname.startswith("Bast_"))] #(rowname.startswith("Ruzicka_MtSinai_") or rowname.startswith("Ruzicka_McLean_"))]
    DF = DF.loc[row_ids_keep]

    #drop columns that are Bast et al, Batiuk res2, Ruzicka Mt Sinai or McLean:
    col_list = DF.columns.tolist()
    cols_drop = [col for col in col_list if col.startswith("Ruzicka_MtSinai_") or col.startswith("Ruzicka_McLean_") or col.startswith("Bast_")]
    DF = DF.drop(columns=cols_drop)
    
    #reorder columns
    col_first = [col for col in DF.columns.tolist() if col.startswith("Batiuk_")]
    col_last = [col for col in DF.columns.tolist()  if col.startswith("Ruzicka_")]
    last_inh_col = np.max([id for id,c in enumerate(col_last) if c.startswith("Ruzicka_In-")])
    col_last = col_last[2:last_inh_col+1]+col_last[0:2]+col_last[last_inh_col+1:len(col_last)]
    DF = DF[col_first + col_last]

    #set values to 0 for Bast vs Bast
    if any(DF.columns.str.startswith("Bast_")):
        row_ids_bast = [rowname for rowname in DF.index.tolist() if rowname.startswith("Bast_")]
        cols_bast = [col for col in DF.columns.tolist()  if col.startswith("Bast_")]
        DF.loc[row_ids_bast,cols_bast] = 0

    return DF

def get_corr_matrix(DF_p,DF_l,opt_corr):
    #correlation matrixes:
    #at least 3 genes overlap to calculate correlation
    if DF_p.index.name != "gene":
        DF_p.set_index("gene",inplace=True)
    DF_p_corr = DF_p.corr(opt_corr,min_periods=3)
    DF_p_corr = make_DF_concise(DF_p_corr)

    if DF_l.index.name != "gene":
        DF_l.set_index("gene",inplace=True)
    DF_l_corr = DF_l.corr(opt_corr,min_periods=3)
    DF_l_corr = make_DF_concise(DF_l_corr)

    return DF_p_corr, DF_l_corr

def prepare_DF_for_upset_plot(DF_p,alpha):
    DF = DF_p.copy()

    col_list = DF.columns.tolist()
    #drop some Ruzicka columns
    cols_drop = [col for col in col_list if col.startswith("Ruzicka_MtSinai_") or col.startswith("Ruzicka_McLean_")]
    DF = DF.drop(columns=cols_drop+["gene"])

    #sign p-val for cutoff alpha?
    DF = DF<=alpha

    #group columns
    col_names = DF.columns.tolist()
    Bast_cols = [c for c in col_names if c.startswith("Bast")]
    Bast_cols_Inh = [c for c in Bast_cols if "Inhibitory" in c]
    Bast_cols_Exc = [c for c in Bast_cols if "Excitatory" in c]
    Bast_cols_NN = [c for c in Bast_cols if c not in Bast_cols_Inh+Bast_cols_Exc]

    Batiuk_cols = [c for c in col_names if c.startswith("Batiuk")]
    Batiuk_cols_Inh = [c for c in Batiuk_cols if "LAMP5" in c or "SST" in c or "VIP" in c or "PVALB" in c or "ID2" in c]
    Batiuk_cols_Exc = [c for c in Batiuk_cols if (("_L2" in c or "_L3" in c or "_L4" in c or "_L5" in c or "_L6" in c) and c not in Batiuk_cols_Inh)]
    #no non neuronal cell types in Batiuk

    Ruzicka_cols = [c for c in col_names if c.startswith("Ruzicka")]
    Ruzicka_cols_Inh = [c for c in Ruzicka_cols if "_In-" in c]
    Ruzicka_cols_Exc = [c for c in Ruzicka_cols if "_Ex-" in c]
    Ruzicka_cols_NN = [c for c in Ruzicka_cols if c not in Ruzicka_cols_Inh+Ruzicka_cols_Exc]

    #add cell class columns to DF
    DF["Bast_Inh"] = DF[Bast_cols_Inh].sum(axis=1)
    DF["Bast_Exc"] = DF[Bast_cols_Exc].sum(axis=1)
    DF["Bast_NN"] = DF[Bast_cols_NN].sum(axis=1)

    DF["Batiuk_Inh"] = DF[Batiuk_cols_Inh].sum(axis=1)
    DF["Batiuk_Exc"] = DF[Batiuk_cols_Exc].sum(axis=1)

    DF["Ruzicka_Inh"] = DF[Ruzicka_cols_Inh].sum(axis=1)
    DF["Ruzicka_Exc"] = DF[Ruzicka_cols_Exc].sum(axis=1)
    DF["Ruzicka_NN"] = DF[Ruzicka_cols_NN].sum(axis=1)

    #drop the other columns:
    DF = DF.drop(columns=Bast_cols + Batiuk_cols + Ruzicka_cols)

    return DF

def plot_upset_gene_overlap_per_class(DF_p,alpha,results_path,min_subset_size,opt_save):
#def plot_number_sign_genes_overlapping_upset(DF,opt_save,path_results,alpha_val,input_data,min_subset_size):
    DF = prepare_DF_for_upset_plot(DF_p,alpha)
    params = {'font.size':10, 'axes.labelsize':10, 'xtick.labelsize':10, 'ytick.labelsize':10}
    with plt.rc_context(params):
        plot(from_indicators(DF.columns.tolist(),data=DF>0), show_counts=True, sort_by = 'cardinality', sort_categories_by='cardinality', min_subset_size=min_subset_size, facecolor="darkgrey")
                
    if opt_save:
        isExist = os.path.exists(results_path)
        if not isExist:
            os.makedirs(results_path)
        plt.savefig(results_path+'upset_DEGs_overlap_within_Cell_classes_'+str(alpha).replace('.','_')+'.pdf', bbox_inches='tight')
    plt.close()
    plt.clf() 