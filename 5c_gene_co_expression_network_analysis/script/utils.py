# author: Lisa Bast
# date: 2024-01-12,  13:46:21
# version: 0.0.1
# about: utils for gene co-expression network analysis

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import plot, from_indicators
import plotly.graph_objects as go
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sctriangulate.colors import build_custom_continuous_cmap

def save_modules_pathways(GO_strings,path_modules_raw,CT_list_modules,opt_test,opt_modules):
    for GO_str in GO_strings:
        file_parent_GO_term_list = path_modules_raw + "reduced_terms_rrvgo_all_"+GO_str+"_size_all_GRN_modules_gprofiler.csv"

        df_parent_term_mapping = pd.read_csv(file_parent_GO_term_list,usecols=["go",'parent',"term","parentTerm"])
        df_parent_term_mapping = get_term_category(df_parent_term_mapping,path_modules_raw)
        #if categories need to be annotated manually, to prepare data frame:
        #df_parent_term_mapping.sort_values(by="parentTerm").to_excel(path_modules_raw+'Term_category_mapping_to_do.xlsx')

        #df_parent_term_mapping["count"] = 1
        #category_abundances = df_parent_term_mapping[["term","category","count"]].groupby(["category","term"]).count()
        df_parent_term_mapping["category"].value_counts()
        df_parent_term_mapping.to_csv(path_modules_raw+"df_parent_term_mapping_"+GO_str+".csv")

        bool_first_ct2 = True
        DF_modules = []
        for ct in CT_list_modules:
            path_modules = path_modules_raw + ct + "/norm_basic/"
            if opt_modules == "scz_implicated":
                module_names = get_scz_relevant_module_names(ct,path_modules_raw,opt_test,0.05)
                for module_name in module_names:
                    module_filename = "module_list_rrvgo_"+module_name+"_"+GO_str+"_size_all.csv"
                    DF_modules, bool_first_ct2 = put_module_results_together(module_filename, module_name, DF_modules,GO_str,path_modules,bool_first_ct2,ct)
            else:
                filenames = os.listdir(path_modules)
                module_filenames = [f for f in filenames if f.startswith("module_list_rrvgo_") and GO_str in f]
                for module_filename in module_filenames:
                    module_name = module_filename.split('_')[3]
                    if module_name =="all":
                        continue
                    DF_modules, bool_first_ct2 = put_module_results_together(module_filename, module_name, DF_modules,GO_str,path_modules,bool_first_ct2,ct)

        DF_modules.drop(columns=[c for c in DF_modules.columns if c not in ["go","Celltype_Module"]],inplace=True)
        DF_modules = pd.merge(left = df_parent_term_mapping, right = DF_modules, how="right",on="go")
        #remove nans
        #DF_modules.dropna(inplace=True)
        if opt_modules=="scz_implicated":
            DF_modules.to_csv(path_modules_raw+"DF_sign_modules_"+GO_str+"_all_cell_types_"+opt_test+".csv")
        else:
            DF_modules.to_csv(path_modules_raw+"DF_all_modules_"+GO_str+"_all_cell_types_"+opt_test+".csv")

def get_scz_relevant_module_names(ct,path_GRN_results,opt_test):
    #open file for respective ct and test
    DME = pd.read_csv(path_GRN_results + ct + "/norm_basic/" + "DMEs_"+opt_test+".csv")
    #determine sign modules
    module_names = DME[DME["p_val_adj"]<=0.05]["module"].tolist()

    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    return module_names 


def get_CT_colors_as_df(n_clusters,path_filtered_data):
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

    return DF

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

def get_module_correlations_and_percentage_overlapping_genes(M_df):

    Correlations = np.zeros(shape=(np.shape(M_df)[0],np.shape(M_df)[0]))
    P_overlapping_genes = np.zeros(shape=(np.shape(M_df)[0],np.shape(M_df)[0]))

    M_df.set_index("module_name_short",inplace=True)
    M_df.drop(['celltype_short',"module"],axis=1,inplace=True)
    for row in range(0,np.shape(M_df)[0]):
        module_name_row = M_df.iloc[row].name
        for col in range(0,np.shape(M_df)[0]):
            if row==col:
                Correlations[row,col]  = 1
                P_overlapping_genes[row,col] = 100
            else:
                module_name_col = M_df.iloc[col].name
                DF_tmp = pd.concat([M_df.iloc[row],M_df.iloc[col]],axis=1)
                #drop all zero entries
                DF_tmp.drop(DF_tmp[DF_tmp.sum(axis=1) == 0].index , inplace=True)
                Correlations[row,col] = DF_tmp[module_name_row].corr(DF_tmp[module_name_col])
                P_overlapping_genes[row,col] = 100*(np.sum(DF_tmp.sum(axis=1)==2)/np.min([DF_tmp[module_name_row].sum(),DF_tmp[module_name_col].sum()]))
    DF_C = pd.DataFrame(Correlations, index = M_df.index, columns=M_df.index.tolist())
    DF_P = pd.DataFrame(P_overlapping_genes, index = M_df.index, columns=M_df.index.tolist())
    DF_C[DF_C.isna()]=0
    DF_P[DF_P.isna()]=0

    return DF_C,DF_P

def get_scz_relevant_module_names(ct,path_GRN_results,opt_test,p_val_cutoff):
    #open file for respective ct and test
    DME = pd.read_csv(path_GRN_results + ct + "/norm_basic/" + "DMEs_"+opt_test+".csv")
    #determine sign modules
    module_names = DME[DME["p_val_adj"]<=p_val_cutoff]["module"].tolist()
    return module_names 

def get_DEP_DEG_module_dataframe(DF,path_genelists,path_module_genes,path_DEG_result, DEG_result_filename,path_results_proteomics,filename_proteomics,alpha_val_DEPs,path_results_proteomics_hgnc,hgnc_file,bordeaux_modules):
    DF_p = get_DEP_results(path_results_proteomics,filename_proteomics,1 ,"","", alpha_val_DEPs ,path_results_proteomics_hgnc+hgnc_file)

    #for the bordeaux modules: plot log2FC in proteomics and highlight MT genes and nuclear genes
    DF_BM = DF[DF["module_name"].isin(bordeaux_modules)]

    DF_BM_p = pd.merge(left=DF_BM, right=DF_p, left_on="accession",right_on="Ensembl gene ID",how="left")
    DF_BM_p.dropna(inplace=True) # not all genes in bordeaux modules map to proteins in DEP data set --> 89 unique genes
    DF_BM_p["MT_gene_status"]=DF_BM_p["Gene"].str.startswith("MT-")

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

def get_DIUG_dataframe(DIUG_filename,alpha_cutoff):
    df_DIG = pd.read_table(DIUG_filename,index_col=0)
    df_DIG["bool_DIG_tested"] = True
    df_DIG["bool_DIG"] = False
    df_DIG.loc[df_DIG["FDR"]<=alpha_cutoff,"bool_DIG"] = True
    df_DIG.rename(columns={"associated_gene": "Gene_short"},inplace=True)

    return df_DIG

def put_module_results_together(module_filename, module_name, DF_modules,GO_str,path_modules,bool_first_ct2,ct):
    if os.path.isfile(path_modules+module_filename):
        df_module = pd.read_csv(path_modules+module_filename)
        df_module["Celltype_Module"] = ct+'_'+module_name
        if len(df_module)>0:
            if bool_first_ct2:
                DF_modules = df_module
                bool_first_ct2 = False
            else:
                DF_modules = pd.concat([DF_modules,df_module])
    else:
        print("no_rrvgo_file_for: " + GO_str + '_', ct+'_'+module_filename)
    return DF_modules, bool_first_ct2

def get_CT_colors_as_dict(n_clusters,path_filtered_data,opt_short):
    if path_filtered_data.endswith('cellranger/'):
        df = pd.read_excel(path_filtered_data+'Cell_type_colors.xlsx',engine='openpyxl')
        if n_clusters==3:
            str_ = 'class'
        else:
            str_ = 'type'
        CT_col = 'Cell '+str_+' ('+str(n_clusters)+')'
        Col_col = 'Cell '+str_+' ('+str(n_clusters)+') color'
        if opt_short:
            df = get_short_CT_names(df,CT_col)
            CT_col_dict = dict(df[[CT_col+"_short",Col_col]].values)
        else:
            CT_col_dict = dict(df[[CT_col,Col_col]].values)
    else:
        print('This is not implemented yet.')
        CT_col_dict = []
    return CT_col_dict 

def get_genes_for_modul_selection(modules,CT_list,path_module_genes):
    bool_first=True
    for m in modules:
        #determine cell type
        ct = [ct for ct in CT_list if m.startswith(ct+"_module_")]
        #get genes for module
        DF_module_genes = pd.read_csv(path_module_genes+"module_gene_mapping_info_"+ct[0]+".csv")
        #filter for relevant modules
        DF_module_genes_sel = DF_module_genes[DF_module_genes["module_name"].isin(modules)]
        if bool_first:
            DF = DF_module_genes_sel.copy()
            bool_first=False
        else:
            DF = pd.concat([DF,DF_module_genes_sel])

    DF["module_name"] = DF["module_name"].str.replace("_module_","_")

    return DF

def get_bordeaux_module_genes_with_pathways_per_category(path_GO_results_BM,path_module_results,path_genes_matrix):
    BM_genes =pd.read_csv(path_GO_results_BM+"df.gene_go_bordeauxModules.tsv",delimiter="\t")
    Module_GOs = pd.read_csv(path_module_results +"DF_sign_modules_BP_all_cell_types_LR.csv")
    DF_BM_genes = pd.merge(right=BM_genes, left =Module_GOs,how = "right",on=["Celltype_Module","go"])
    #genematrix: get gene names:
    GM = pd.read_csv(path_genes_matrix + "geneMatrix_v2.tsv",delimiter="\t")
    DF_BM_genes = pd.merge(right=DF_BM_genes, left = GM[["ensgid","gene_name"]], right_on="ENSGID",left_on = "ensgid",how="right")
    DF_BM_genes.drop(columns=["Unnamed: 0"],inplace=True)
    #filter for genes in bordeaux modules:
    genes_in_bordeaux_modules = pd.read_csv(path_genes_matrix +"magenta_module_genes_gene_ids.csv",header=None,names=["gene_name"])
    DF_BM_genes = pd.merge(right=DF_BM_genes, left = genes_in_bordeaux_modules, on="gene_name",how="inner")

    return DF_BM_genes

def get_hub_gene_DF(CT_list, path_modules_raw,modules,mapping):
    bool_first=True
    #for ct in CT_list:
    for m in modules:
        ct = m.split('_module_')[0]
        module = m.split('_module_')[1]
        df_hub = pd.read_csv(path_modules_raw+ct+"/norm_basic/"+"hub_df_all.csv")
        df_hub = df_hub[df_hub["module"]==module]
        df_hub.drop(columns="Unnamed: 0",inplace=True)
        df_hub["celltype"] = ct
        df_hub["module_name"] = ct+'_'+df_hub["module"]
            
        if bool_first:
            DF_hub=df_hub.copy()
            bool_first=False
        else:
            DF_hub=pd.concat([DF_hub,df_hub])
    
    #filter for relevant modules:
    #DF_hub = DF_hub[DF_hub["module_name"].isin(modules)]
    #add ensgid:
    DF_hub = pd.merge(left=DF_hub,right=mapping[["ensgid","gene_name"]],on="gene_name",how="left")

    return DF_hub

def get_module_selection(opt_module_selection,CT_list,path_module_genes,opt_test):
    if opt_module_selection == "bordeaux":
        modules = ["Astrocytes_module_yellow",
                            "Excitatory_Layer_5-6_IT_neurons_I_module_green",
                            "Oligodendrocytes_module_red",
                            "Inhibitory_LAMP5_neurons_module_red",
                            "Inhibitory_VIP_neurons_module_pink",
                            "Excitatory_Layer_2-3_IT_neurons_I_module_red",
                            "Excitatory_Layer_2-3_IT_neurons_II_module_black", 
                            "Inhibitory_PVALB_neurons_module_magenta"]
    elif opt_module_selection == "midnight":
        modules = ["Excitatory_Layer_3-6_IT_neurons_module_black", 
                            "Inhibitory_LAMP5_neurons_module_yellow", 
                            "Oligodendrocytes_module_yellow", 
                            "Excitatory_Layer_2-3_IT_neurons_I_module_blue",
                            "Excitatory_Layer_2-3_IT_neurons_II_module_blue"]
    else:
        modules = []
        for ct in CT_list:
            if opt_module_selection == "significant_scz":
                module_color = get_scz_relevant_module_names(ct,path_module_genes,opt_test)
            elif opt_module_selection == "all":
                module_color = get_module_names(ct,path_module_genes,opt_test)
            modules_ct = [ct+'_'+m for m in module_color]
            modules = modules + modules_ct
    return modules

def get_module_names(ct,path_GRN_results,opt_test):
    #open file for respective ct and test
    DME = pd.read_csv(path_GRN_results + ct + "/norm_basic/" + "DMEs_"+opt_test+".csv")
    #determine sign modules
    module_names = DME["module"].tolist()

    #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
    return module_names 

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

def get_term_category(df_term_mapping,path_modules_raw):

    df_term_category_mapping = pd.read_excel(path_modules_raw+"Term_category_mapping.xlsx",sheet_name="rrvgo_child_parent_terms")
    df_term_mapping = pd.merge(left=df_term_mapping,right=df_term_category_mapping,on="term",how="left")

    return df_term_mapping

def plot_sankey_pathways_2_levels(df,GO_str, CT_class,DF_colors,opt_no_pw_color):
    nodes = np.unique(df[["source","target"]], axis=None)
    nodes = pd.Series(index=nodes, data=range(len(nodes)))
    
    #colors:
    colors = nodes.index.tolist()
    #colors = [c.replace("_"," ") for c in colors]
    for i,c in enumerate(colors):
        if opt_no_pw_color:
            if c.rsplit(' ', 1)[0] in DF_colors["celltype_short"].tolist():
                colors[i] = DF_colors[DF_colors["celltype_short"]==c.rsplit(' ', 1)[0]]["color"].tolist()[0]
            else:
                colors[i] = "#808080"
        else:
            if "RNA processing" in c or "RNA splicing" in c or "mRNA processing" in c or "interstrand cross-link repair" in c or "mRNA metabolic process" in c or "transcription" in c or "gene expression" in c or "spliceo" in c or "nucleus" in c or "chromosome organi" in c or "nucleobase-containing" in c or "nucleotide" in c or "nucleoside" in c or "RNA metab" in c:
                colors[i] = "#7F00FF"
            elif c.rsplit(' ', 1)[0] in DF_colors["celltype_short"].tolist():
                colors[i] = DF_colors[DF_colors["celltype_short"]==c.rsplit(' ', 1)[0]]["color"].tolist()[0]
            elif "mitochon" in c or "Mitoch" in c or "ATP" in c or "biosynthetic process" in c or "aerobic respiration" in c or "oxidative phosphory" in c or "electron transport" in c or "energy" in c:
                colors[i] = "#9E0B0B"
            elif "protein folding" in c or "protein refolding" in c or "protein modification process" in c or "protein maturation" in c  or "protein phosphorylation" in c or "protein secretion" in c  or "translation" in c:
                colors[i] = "#0B469E" 
            # elif "GABA" in c: 
            #     df_parent_term_mapping.at[index,"category"] = "GABA cell signaling"
            # elif "AMPA" in c or "glutamate" in c:
            #     df_parent_term_mapping.at[index,"category"] = "Glutamate cell signaling"
            # elif "calcium" in c or "catecholamine" in c:
            #     df_parent_term_mapping.at[index,"category"] = "calcium cell signaling"
            elif "GABA" in c or "AMPA" in c or "glutamate" in c or "calcium" in c or "catecholamine" in c or "synap" in c or "cell signaling" in c or "receptor" in c or "neurotransmitter" in c or "vesicle" in c or " ion" in c or c=="ion transport" or "ion homeostasis" in c or "cation" in c or "cell communication" in c or "hormone" in c:
                colors[i] = "#EC9600"
            elif "transmembrane trans" in c or "amino acid transport" in c or "chloride transport" in c or "hormone transport" in c or "import into cell" in c or "export from cell" in c or "cell adhesion" in c or "channel complex" in c or "signal transduction" in c or "cell junction" in c or "membrane" in c:
                colors[i] = "#EC9600"
            elif "inflammatory response" in c or "T cell " in c or "cytokine" in c or "immune system" in c:
                colors[i] = "#808080"
            elif "memory" in c or "learning" in c:
                colors[i] = "#808080"
            elif "neuron projection" in c or "dendrite" in c or "axon" in c or "myelin" in c:
                colors[i] = "#808080"
            elif "muscle" in c:
                colors[i] = "#808080"
            elif "behavior" in c or "locomotion" in c:
                colors[i] = "#808080"
            elif "apopto" in c or "differentiation" in c or "homeostatic process" in c or "homeostasis" in c or "generation of neurons" in c or "proliferation" in c or "cell cycle" in c or "migration" in c or "tissue remodeling" in c or "cell motility" in c or "cell death" in c or "cell fate " in c or "gliogenesis" in c or "heart growth" in c or "cell activation" in c or "maturation" in c or "organ growth" in c or "neurogenesis" in c or " stem cell" in c or "axonogenesis" in c or "angiogenesis" in c or "cell killing" in c or "heart growth" in c or "development" in c or "hemopoiesis" in c or "morphogenesis" in c or "cell growth" in c or "mitotic" in c or "mitosis" in c:
                colors[i] = "#009999"
            else:
                colors[i] = "#808080"

        
    #add cell type colors first:
    #colors = colors.replace(colors.startswith("Excitatory_Layer_2-3_IT_neurons_II"))
    fig = go.Figure(data=[go.Sankey(
        node = dict(
            pad = 15,
            thickness = 20,
            line = dict(color = "black", width = 0.5),
            label = nodes.index,
            color = colors
        ),
            link = dict(
            source = nodes.loc[df["source"]],
            target = nodes.loc[df["target"]],
            value = df["value"]
        ))]
    )

    fig.update_layout(title = "Pathway Parent Terms GO:"+GO_str+" "+CT_class+" modules", font_size=10)
    fig.show()
    #plot(fig, auto_open=True)

def plot_module_pathway_heatmap_and_clustermap(DF_modules,GO_str,  df_mapping, path_results, subfolder,opt_test,opt_dbscan,opt_order_by_cluster,opt_order_by_ct):

    modules_ordered = pd.read_csv(path_results+"modules_ordered_"+opt_test+".csv")
    modules_ordered["module"] = modules_ordered["module"].str.replace("L ","L") 
    modules_ordered.drop(columns="Unnamed: 0", inplace=True)
    path_results = path_results + subfolder
    if opt_dbscan==True:
        variab = "cluster"
        fig1_w = 8
        fig1_h = 4
        fig2_w = 8
        fig2_h = 6
        color_map = "Blues"
    else:
        variab = "category"
        fig1_w = 11
        fig1_h = 6
        color_map = "Reds"
    DF_modules['number pathways']=1


    clustermap_data = pd.pivot_table(DF_modules, values="number pathways", index=variab, columns = 'Celltype_Module_short',aggfunc="sum")
    #check if columns are missing and if so add with zero:
    for m in modules_ordered['module'].tolist():
        if m not in clustermap_data.columns:
            clustermap_data[m]=0
    clustermap_data = clustermap_data[modules_ordered['module'].tolist()]
    clustermap_data.fillna(0, inplace=True)
    #heatmap_data_normalize_by_number_pws_in_category = heatmap_data.div(heatmap_data.sum(), 1)
    #heatmap_data_normalize_by_number_pws_in_category_and_cluster = heatmap_data_normalize_by_number_pws_in_category.div(heatmap_data_normalize_by_number_pws_in_category.sum(1), 0)
    clustermap_data_normalize_by_number_pws_in_category = clustermap_data.div(df_mapping[variab].value_counts(), 0)
    #sns.heatmap(heatmap_data_normalize_by_number_pws_in_category, cmap="Blues")
    #heatmap_data_normalize_by_number_pws_of_ct = heatmap_data.div(heatmap_data.sum(0), 1)
    #sns.heatmap(heatmap_data_normalize_by_number_pws_of_ct, cmap="Blues")
    clustermap_data_normalize_by_number_pws_in_category_and_of_ct = clustermap_data_normalize_by_number_pws_in_category.div(clustermap_data_normalize_by_number_pws_in_category.sum(0), 1)
    clustermap_data_normalize_by_number_pws_in_category_and_of_ct.fillna(0, inplace=True)
    #fig, axes  = plt.subplots(1, 4, sharey='row', gridspec_kw={'width_ratios':[1,6,5,1]})
    #fig2, axes2 = plt.subplots()
    if opt_dbscan==True:
        cg = sns.clustermap(clustermap_data_normalize_by_number_pws_in_category_and_of_ct, cmap=color_map,xticklabels=True, yticklabels=True,figsize=(fig1_w,fig1_h), col_cluster = True, row_cluster=True)#,ax=axes[1],cbar_ax=axes[0])
    else:
        categories_ordered = ["Synapse",
                              "Cell adhesion",
                              "Signal transduction",
                              "Other Development",
                              "Nervous system development",
                              "Kinase/Peptidase",
                              "Ion transport",
                              "Energy (ATP)",
                              "Protein modification",
                              "Protein assembly",
                              "Extracellular matrix",
                              "Filament organization",
                              "Gene expression",
                              "Cellular process",
                              "Cell-cell signalling",
                              "Calcium signalling",
                              "Other transport",
                              "Cellular response to stimuli",
                              "Female reproduction",
                              "Intracellular signal transduction",
                              "Angiogenesis/Blood",
                              "Apoptosis",
                              "Inflammation"]

        clustermap_data_normalize_by_number_pws_in_category_and_of_ct = clustermap_data_normalize_by_number_pws_in_category_and_of_ct.reindex(categories_ordered)
        cg = sns.clustermap(clustermap_data_normalize_by_number_pws_in_category_and_of_ct, cmap=color_map,xticklabels=True, yticklabels=True, col_cluster = True, row_cluster=False, figsize=(fig1_w,fig1_h))#,ax=axes[1],cbar_ax=axes[0])
    plt.savefig(path_results+'Clustermap_module_pathways_per_'+variab+'_'+GO_str+'_'+opt_test+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

    if opt_dbscan==True:
        if opt_order_by_cluster:
            #get order of clustermap y-values:
            cluster_order = clustermap_data_normalize_by_number_pws_in_category_and_of_ct.index[cg.dendrogram_row.reordered_ind].tolist()
        heatmap_data = pd.pivot_table(DF_modules, values="number pathways", index="cluster", columns = 'category',aggfunc="sum")
        heatmap_data.fillna(0, inplace=True)
        heatmap_data_normalize_by_number_pws_in_cluster = heatmap_data.div(heatmap_data.sum(1), 0)
        heatmap_data_normalize_by_number_pws_in_category_and_cluster = heatmap_data_normalize_by_number_pws_in_cluster.div(heatmap_data_normalize_by_number_pws_in_cluster.sum(0), 1)
        if opt_order_by_cluster:
            heatmap_data_normalize_by_number_pws_in_category_and_cluster = heatmap_data_normalize_by_number_pws_in_category_and_cluster.reindex(cluster_order,level="cluster")
        #
        #fig, ax = plt.subplots(figsize=(8,4))  
        #sns.heatmap(heatmap_data_normalize_by_number_pws_in_category_and_cluster, cmap="Blues", square=True)#,ax=axes[2],cbar_ax=axes[3])
        sns.clustermap(heatmap_data_normalize_by_number_pws_in_category_and_cluster, cmap="Reds",xticklabels=True, yticklabels=True, col_cluster = True, row_cluster=False, figsize=(fig2_w,fig2_h))#,ax=axes[2],cbar_ax=axes[3])
        plt.savefig(path_results+'Heatmap_cluster_category_'+GO_str+'_'+opt_test+'.pdf',bbox_inches='tight')
        plt.close()
        plt.clf()
    else:
        # plot number of co-occurrence of biol. functions in clusters
        var_of_interest = "Celltype_Module_short"#"cluster"
        df = DF_modules[["category",var_of_interest]].copy()
        df = df.drop_duplicates().sort_values(var_of_interest)
        categories = df["category"].unique().tolist()
        M = np.zeros(shape=(len(categories),len(categories)))
        for cluster in df[var_of_interest].unique().tolist():
            cat_sel = df[df[var_of_interest]==cluster]["category"].tolist()
            idx = np.argwhere(np.isin(categories, cat_sel)).ravel()
            for j in idx:   
                for i in idx:
                    M[j,i] = M[j,i]+1 

        #normalize by number of occurrences of the pathway:
        N = df["category"].value_counts()
        DF_M = pd.DataFrame(data=M,
                            index=categories,
                            columns=categories)
        for m in range(0,len(categories)):#columns
            DF_M_norm = DF_M.copy()
            DF_M_norm[categories[m]] = DF_M_norm[categories[m]]/N[categories]
        for n in range(0,len(categories)):#rows
            DF_M_norm.iloc[n] = DF_M_norm.iloc[n]/N[categories]    

        #N[categories[i]]
        sns.heatmap(DF_M_norm,cmap="Blues")
        plt.savefig(path_results+'Heatmap_category_cooccurrence_in_'+var_of_interest+'_normalized_'+GO_str+'_'+opt_test+'.pdf',bbox_inches='tight')
        plt.close()
        plt.clf()

        cm = sns.clustermap(DF_M,cmap="Blues")
        plt.close()
        category_order = DF_M.index[cm.dendrogram_row.reordered_ind].tolist()
        DF_M_reordered = DF_M.reindex(category_order)
        n_modules = len(df[var_of_interest].unique().tolist())
        mask = np.tril(np.ones_like(DF_M_reordered[category_order], dtype=bool))
        mask[np.diag_indices_from(mask)] = False

        sns.heatmap(DF_M_reordered[category_order],cmap="Blues",mask=mask,vmin=1, vmax = len(category_order),cbar_kws={'ticks': [t for t in range(1,n_modules)]},square=True)
        plt.savefig(path_results+'Heatmap_category_cooccurrence_in_'+var_of_interest+'_'+GO_str+'_'+opt_test+'.pdf',bbox_inches='tight')
        plt.close()
        plt.clf()

        #sort processes according to hm labels
        if opt_dbscan==False:
            pw_categories_per_module = clustermap_data_normalize_by_number_pws_in_category_and_of_ct.reindex(category_order, axis=0)
            if opt_order_by_ct==True:
                columns_ordered_by_ct = sorted(pw_categories_per_module.columns.tolist())
                pw_categories_per_module = pw_categories_per_module[columns_ordered_by_ct]

            sns.heatmap(pw_categories_per_module, cmap=color_map,xticklabels=True, yticklabels=True, square=True)#,ax=axes[1],cbar_ax=axes[0])
            plt.savefig(path_results+'Heatmap_module_pathways_per_'+variab+'_'+GO_str+'_'+opt_test+'_sorted.pdf',bbox_inches='tight')
            plt.close()
            plt.clf()

def plot_number_sign_pathways_as_upset_modules(DF_modules,opt_test,GO_str,path_results,min_ss,opt_save):
    params = {'font.size':10, 'axes.labelsize':10, 'xtick.labelsize':10, 'ytick.labelsize':10}
    color='grey'
    DF_modules["value"] = 1
    DF_bool_red = pd.pivot_table(data = DF_modules, index=["go"], columns = "Celltype_Module_short", values = "value", fill_value=0).astype("bool")
    modules = DF_modules["Celltype_Module_short"].unique().tolist()
    
    #min subset size in upset plot
    #at least one CT needs to have more than min_ss
    # at least one pathway needs to be implicated in more than one cell type
    if np.max(DF_bool_red.sum(axis=0))>min_ss and np.any(DF_bool_red.sum(axis=1)>1):
        #rename columns (add spaces)
        column_names = DF_bool_red.columns.tolist()
        for col in column_names:
            eval("DF_bool_red.rename(columns={\'"+col+"\':\'               "+col+"\'},inplace=True)")
        # update CTs
        CTs = DF_bool_red.columns.tolist()
        if len(CTs)>np.max([2,min_ss]):
            with plt.rc_context(params):
                plot(from_indicators(CTs,data=DF_bool_red), show_counts=True, sort_by = 'cardinality', sort_categories_by='cardinality', min_subset_size=min_ss, facecolor=color)# 
            if opt_save:
                path=path_results+'figures/upset/'
                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                plt.savefig(path+'upset_significant_pathways_'+GO_str+'_'+opt_test+'.pdf', bbox_inches='tight')
            plt.close()
            plt.clf() 

def plot_ranked_log2FC(DEG_MT_sorted,fig_title,path_results,hue_var,n_cluster,path_colors):
    if hue_var == "regulation":
        color_dict = dict(down="teal",up="darkred",unchanged="grey")
    else:
        color_dict = get_CT_colors_as_dict(n_cluster,path_colors,False)
    p = sns.scatterplot(data=DEG_MT_sorted,x="x1", y="log2FoldChange",hue = hue_var,palette=color_dict,s=13)
    p.axhline(y=0,color="black") 
    p.set_title(fig_title)
    p = remove_frame_keep_axes(p)
    plt.savefig(path_results+"Ranked_log2FC_"+fig_title.replace(" ","_")+"_colored_by_"+hue_var+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf()

def plot_proteomics_vs_transcriptomics(DF,opt_filtering,path_results):
    for s in ["log2FoldChange","signed adjusted p-value"]:
        fig,ax = plt.subplots(1,1,figsize=(5,5))
        if opt_filtering=="MT_genes_filtered":
            sns.scatterplot(data=DF[DF["MT_gene_status"]==True],y=s+"_DEP",x=s+"_DEG",s=5,ax=ax)
        elif opt_filtering=="nuclear_genes_filtered":
            sns.scatterplot(data=DF[DF["MT_gene_status"]==False],y=s+"_DEP",x=s+"_DEG",s=5,ax=ax)
        else:
            sns.scatterplot(data=DF,y=s+"_DEP",x=s+"_DEG",s=5,ax=ax)
        ax.axhline(y=0,color="lightgrey")
        ax.axvline(x=0,color="lightgrey")
        fig.savefig(path_results+s+"_proteomics_transcriptomics_Bordeaux_Module_genes_"+opt_filtering+".pdf",bbox_inches='tight')
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

def plot_number_genes_per_module(M_df,n_cl,path_module_genes,path_filtered_loom_data,opt_test,opt_filter_genes):
    pd.options.display.float_format = '{:.2f}'.format
    colors = get_CT_colors_as_df(n_cl,path_filtered_loom_data)
    colors = get_short_CT_names(colors,"celltype")
    colors["celltype_short"] = colors["celltype_short"].str.replace("_"," ")
    M_df.reset_index(inplace=True)
    M_df[["celltype_short","module"]] = M_df["module_name_short"].str.rsplit(pat=" ",n=1,expand=True)
    M_df = pd.merge(left = M_df, right=colors,on="celltype_short")
    #sort cts:
    #M_df["module_name_short"] =
    M_df.set_index(keys="module_name_short",drop=True, inplace=True)
    M_df.drop(columns=["celltype","celltype_short","module"],inplace=True)

    M_df_ordered = M_df.loc[M_df.index[M_df.index.str.startswith("Exc")].tolist() + M_df.index[M_df.index.str.startswith("Inh")].tolist() + M_df.index[np.logical_and(M_df.index.str.startswith("Exc")==False,M_df.index.str.startswith("Inh")==False)].tolist()]
    pl = M_df_ordered.loc[:,M_df_ordered.columns != "color"].sum(axis=1).plot(kind="barh",color=M_df_ordered["color"],xlabel = "Number of Genes")
    pl.ticklabel_format(axis="x", style="plain")
    pl = remove_frame_keep_axes(pl)
    plt.savefig(path_module_genes+'Number_genes_per_module'+opt_test+'_'+opt_filter_genes+'.pdf',bbox_inches='tight')
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

def plot_number_nuclear_mito_and_mito_genes_per_bordeaux_module(DF,path_results):
    DF = DF[["gene_name","Celltype_Module_short"]].drop_duplicates()
    DF = DF.reset_index()
    DF["type"] = "nuclear"
    DF.loc[DF["gene_name"].str.startswith("MT-"),"type"] = "mitochondrial"
    DF.drop(columns="index").to_csv(path_results+"nuclear_and_mitochondrial_genes_per_bordeaux_module.csv",header=True)
    DF["count"] = 1
    DF["Celltype_Module_short"] = DF["Celltype_Module_short"].str.replace("_"," ")
    fig, ax = plt.subplots()
    DF_p = pd.pivot_table(data =DF,values="count",index = "Celltype_Module_short", columns = ["type"],aggfunc="sum")
    DF_p.plot(kind="bar",color=["darkgrey","grey"],ax=ax)
    ax = remove_frame_keep_axes(ax)
    plt.savefig(path_results+"Number_of_nuclear_and_mitchondrial_genes.pdf",bbox_inches='tight')
    plt.close()
    plt.clf()
    
def plot_volcano_all_modules(R_DME,n_cl,path_data,path_results,opt_test,opt_save):
    plt = loadPltSettings(12,10)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(7,4))
    fig.tight_layout()
    color_dict = get_CT_colors_as_dict(n_cl,path_data,True)
    R_DME = get_short_CT_names(R_DME,"level")
    sns.scatterplot(data=R_DME,x="avg_log2FC", y="-log10(adj p-val)",hue="level_short",palette=color_dict,ax=ax) 
    #add line for sign th alpha
    ax.axhline(np.log10(0.05)*(-1),c="grey",linewidth=1,linestyle="--")
    ax.axvline(0.0,c="grey",linewidth=1,linestyle="--")
    for ID, yval in enumerate(R_DME["-log10(adj p-val)"]):
        # annotate modules with largest effectsize and modules with smallest p-val
        if np.abs(yval)>35 or (np.abs(R_DME.iloc[ID].avg_log2FC) > 0.5 and np.abs(yval)>np.log10(0.05)*(-1)):
            ax.annotate(R_DME.iloc[ID].module,(R_DME["avg_log2FC"].iloc[ID],yval),fontsize=10)
    #sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax = remove_frame_keep_axes(ax)
    ax.ticklabel_format(axis="x", style="plain")
    ax.ticklabel_format(axis="y", style="plain")
    if opt_save:
        path=path_results+'figures/volcano_DMEs/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        plt.savefig(path+'Volcano_pot_all_DMEs_'+opt_test+'.pdf',bbox_inches='tight')
    plt.close()
    plt.clf() 

def plot_GOs_as_horizontal_barplot(DF,add_str,path_results,variable_term,variable_val,GO_str,n_top_pws):
    #sort DF
    DF = DF.sort_values(by=variable_val,ascending=True)
    if n_top_pws!=999:
        DF = DF.head(n_top_pws)
        add_str_top_pws = str(n_top_pws)
    else:
        add_str_top_pws = ""
    if GO_str=="GO BP" and add_str == "_bordeaux_modules_1_overlapping_genes":
        fsize=(5,13)
    elif GO_str=="SYNGO" and add_str == "_bordeaux_modules_2_overlapping_genes":
        fsize=(5,6)
    elif GO_str=="GO BP" and add_str == "_genes_Excitatory_Layer_2-3_IT_neurons_II_red":
        fsize=(3,5)
    elif GO_str=="SYNGO" and add_str == "_genes_Excitatory_Layer_2-3_IT_neurons_II_red":
        fsize=(2,3)
    else: 
        fsize=(5,10)
    fig,ax = plt.subplots(1,1,figsize=fsize)
    DF.plot.barh(y=variable_val,x=variable_term,ax=ax)
    #ax.invert_yaxis()
    if variable_val=="negative_log10_of_adjusted_p_value" or variable_val=="negative log10 FDR corrected p-value":
        ax.set_xlabel("-log10(adj p-value)")
    ax.set_ylabel(GO_str+" term")
    ax=remove_frame_keep_axes(ax)
    ax.get_legend().remove()
    fig.savefig(path_results+"sign_"+GO_str.replace(" ","_")+add_str_top_pws+"_pws_"+add_str+".pdf",bbox_inches='tight')
    plt.close()
    plt.clf()

def plot_sign_syngo_pws(DF_S,path_results):
    sns.heatmap(data = DF_S,square=True, vmin = 0, vmax=0.06, cmap="Greys_r")#,"GSEA 'gene cluster' FDR corrected p-value")
    plt.savefig(path_results+"sign_SYNGO_pws.pdf",bbox_inches='tight')
    plt.close()
    plt.clf()

#plot connectivity distribution per module, colored by cell type
def plot_connectivity_distribution(DF_hub,opt_module_selection,path_results,mode):
    #create results path if non-existent
    if not os.path.exists(path_results):
        os.makedirs(path_results) 

    if opt_module_selection=="significant_scz":
        y = 15
    elif opt_module_selection=="all":
        y = 30
    else:
        y = 5
    #beeswarm plot
    if mode=="DEGs":
        fig,ax = plt.subplots(2,1,figsize=(y,5),sharex=True)

        DF_hub.sort_values(by="gene_type",inplace=True)
        pal = {"lncRNA":"#F08B18","RBP coding gene" :"#00B3FF", "other protein coding gene":"#1D456C", "MT gene":"#CC0000", "nuclear mito gene": "#329C6E","other":"lightgrey"}#[,"#1D456C","#F08B18","#920F0F"]#colors={"lncRNA":"green","other protein coding gene":"blue","RBP coding gene":"purple","MT gene":"red"}
        sns.swarmplot(data=DF_hub, x="module_name_short",y="kME",hue="gene_type",ax=ax[0],size=3,palette=pal)
        ax[0].tick_params(axis="x",rotation=90)
        ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        DF_hub.sort_values(by="regulation status",inplace=True)
        pal2 = {"not_tested":"lightgrey", "unchanged":"grey", "upregulated":"darkred", "downregulated":"teal"}#["lightgrey",,"teal","darkred"] 
        sns.swarmplot(data=DF_hub, x="module_name_short",y="kME",hue="regulation status",ax=ax[1],size=3,palette=pal2)
        ax[1].tick_params(axis="x",rotation=90)
        ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.savefig(path_results+"connectivitiy_beeswarm_"+opt_module_selection+"_module_genes_"+mode+".pdf",bbox_inches='tight')
        plt.close()
        plt.clf()  

        fig,ax = plt.subplots(1,1,figsize=(3,5))
        p=sns.kdeplot(data=DF_hub,x="kME",hue="module_name_short",fill=True, common_norm=False, alpha=0.2)
        p.set_xlim([DF_hub["kME"].min(),DF_hub["kME"].max()])
        #p.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        sns.move_legend(p, "upper left", bbox_to_anchor=(1, 1))
        p.set_title(opt_module_selection+" modules' connectivity distribution")

        fig.savefig(path_results+"connectivitiy_density_"+opt_module_selection+"_module_genes.pdf",bbox_inches='tight')
        plt.close()
        plt.clf()  

    elif mode=="proteomics":
        limit_cb = np.max(np.abs([DF_hub[DF_hub["log2FoldChange_DEP"].isna()==False]["log2FoldChange_DEP"].min(),
                                  DF_hub[DF_hub["log2FoldChange_DEP"].isna()==False]["log2FoldChange_DEP"].max()]))
        from matplotlib import colors, colorbar
        CM = build_custom_continuous_cmap(np.array(fc.dec_to_rgb(colors.to_rgb('purple'))).tolist(),
                                      np.array(fc.dec_to_rgb(colors.to_rgb('lightgrey'))).tolist(),
                                      np.array(fc.dec_to_rgb(colors.to_rgb('darkorange'))).tolist())
        CM.set_bad("black")
        # Normalize to the range of possible values from df["c"]
        norm = colors.Normalize(vmin=(-1)*limit_cb, vmax=limit_cb)
        # create a color dictionary (value in c : color from colormap) 
        colors = {}
        for cval in DF_hub["log2FoldChange_DEP"]:
            colors.update({cval : CM(norm(cval))})

        fig,ax = plt.subplots(1,1,figsize=(y,5))
        #sns.swarmplot(data=DF_hub, x="module_name_short", y="kME", ax=ax, size=5, color='darkgrey')
        sns.swarmplot(data=DF_hub, x="module_name_short", y="kME", hue="log2FoldChange_DEP", ax=ax, size=5, palette=colors,vmin = (-1)*limit_cb, vmax = limit_cb)
        ax.tick_params(axis="x",rotation=90)
        #ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gca().legend_.remove()
        ## create colorbar ##
        divider = make_axes_locatable(plt.gca())
        ax_cb = divider.new_horizontal(size="5%", pad=0.05)
        fig.add_axes(ax_cb)
        cb1 = colorbar.ColorbarBase(ax_cb, cmap=CM,
                                        norm=norm,
                                        orientation='vertical')
        cb1.set_label("log2FC in DEP analysis")
        plt.show()

        fig.savefig(path_results+"connectivitiy_beeswarm_"+opt_module_selection+"_module_genes_"+mode+".pdf",bbox_inches='tight')
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