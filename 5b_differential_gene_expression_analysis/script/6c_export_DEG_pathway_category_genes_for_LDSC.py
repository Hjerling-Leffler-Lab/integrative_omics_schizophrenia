# author: Lisa Bast
# date: 2024-05-30,  17:06:55
# version: 0.0.1
# about: put together all genes on pathways per rrvgo category for DEG pathways --> basis for testing which categories/ functions are enriched for scz risk (to do: Shuyang)

import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import gseapy as gp

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

main_path = path_code.replace("5b_differential_gene_expression_analysis\\script","")

alpha_vals = [0.1,0.05]
GO_categories = ["CC","BP"]#,"MF"]
reg_strings = ["up","down"]
rrvgo_opt = "scores" #"size"
opt_pw_genes = False#True # weather to use all genes on pathway(True) or filter for DEGs (False)
path_gene_list = main_path + "/5b_differential_gene_expression_analysis/data/Gene_lists/rrvgo_category_DEGs/"
path_GO_DB = main_path + "/5b_differential_gene_expression_analysis/data/Gene_lists/GO_DB/"
path_DEG_PW_results = main_path + "/5b_differential_gene_expression_analysis/output/GSA_analysis/DEGs/v3/"
path_DEG_result = main_path + "/5b_differential_gene_expression_analysis/output/DEGs/design_with_6_RUVs/"
DEG_filename = "T5_DEGs_per_cell_type.csv" #"df_DESeq2_all_genes_all_CTs.csv"

def get_ensgids_for_go_category(GO_ids_tmp,GO_dict):
    genes_list_tmp = []
    for go_id in GO_ids_tmp:
        if go_id in dict.keys(GO_dict):
            genes_list_tmp = genes_list_tmp + GO_dict[go_id]
        else:
            print(go_id + " of category "+cat+ "("+str(len(GO_ids_tmp))+" IDs in total) not in GO data base!")     
    return genes_list_tmp

for alpha in alpha_vals:
    alpha_str = str(alpha).replace("0.","")
    path_deg_pw_result = path_DEG_PW_results + "alpha05/min_overlap_3/6_RUVs/cortex/"
    if alpha!=0.05:
        #load DEG results, get DEGs for alpha threshold
        DF_DEGs = pd.read_csv(path_DEG_result+DEG_filename)
        all_DEG_ensgids_for_alpha = DF_DEGs[DF_DEGs["padj"]<=alpha]["Gene"].str.split("_").str[1].tolist()
    for GO_str in GO_categories:
        #get GO terms for rrvgo category:
        for reg_str in reg_strings:
            file_name_rrvgo_result = "reduced_terms_rrvgo_all_"+GO_str+"_"+rrvgo_opt+"_all.csv"
            df = pd.read_csv(path_deg_pw_result+file_name_rrvgo_result)
            df.drop(columns=["cluster","parent","parentSimScore","score","size"],inplace=True)
            
            #get pw result with ensgids:
            file_name_pw_result = "DF_all_"+reg_str+".csv"
            df2_tmp = pd.read_csv(path_deg_pw_result+file_name_pw_result)
            #keep only sign pathways:
            df2_tmp = df2_tmp[df2_tmp["P.fdr.group"].str.replace(",",".").astype("float")<0.05]
            #keep only current GO term:
            df2_tmp = df2_tmp[df2_tmp["group"].str.endswith(GO_str)]

            #keep only relevnat columns:
            df2_tmp[["GO_ID", "pathway_name"]]=df2_tmp["geneset"].str.split(" ",expand=True)
            df2_tmp = df2_tmp[["genes.in.geneset.tested","GO_ID","geneset","P.fdr.group"]]

            #combine data frames of up and down:
            if reg_str!=reg_strings[0]:
                df2 = pd.concat([df2,df2_tmp])
            else:
                df2 = df2_tmp
        #split column into two geneset:
        df2[["ID","term"]]=df2["geneset"].str.split(" ",n=1,expand=True)
        if opt_pw_genes or alpha!=0.05:
            #load database:
            GO_dict = gp.read_gmt(path_GO_DB+"hsapiens.GO_"+GO_str+".ENSG.gmt")
        #for each rrvgo category:
        rrvgo_categories = df["parentTerm"].unique().tolist()
        for cat in rrvgo_categories:
            print("--------------"+cat+"-------------------")
            print(str(len(rrvgo_categories)) + " rrvgo categories in " + GO_str)
            #unique list of GO terms per rrvgo category
            GO_ids_tmp = df[df["parentTerm"]==cat]["go"].tolist()
            #find relevant GO IDs:
            if opt_pw_genes:
                sub_dir = "all_GO_pw_genes/"
                #find relevant GO terms of this rrvgo category
                genes_list_tmp = get_ensgids_for_go_category(GO_ids_tmp,GO_dict)
            else:
                sub_dir = "filtered_for_DEGs_alpha"+alpha_str+"/"
                if alpha==0.05:
                    df_genes = df2[df2["ID"].isin(GO_ids_tmp)]
                    #find DEGs mapping to those GO terms
                    genes_list_tmp = df_genes["genes.in.geneset.tested"].dropna().str.split("_").tolist()
                    #split up ensgids:
                    genes_list_tmp = sum(genes_list_tmp,[])
                    genes_list_tmp = [g for g in genes_list_tmp if g!="NA"]
                else:
                    #To do: filter DEG list for genes on pathway
                    genes_list_tmp = get_ensgids_for_go_category(GO_ids_tmp,GO_dict)
                    print(str(len(genes_list_tmp)))
                    genes_list_tmp = [g for g in genes_list_tmp if g in all_DEG_ensgids_for_alpha]
                    print(str(len(genes_list_tmp)))
            #create path
            os.makedirs(path_gene_list+"alpha05/"+sub_dir, exist_ok=True)
            #save a file for each rrvgo category
            cat = cat.replace("/"," or ")
            np.savetxt(path_gene_list+"alpha05/"+sub_dir+"GO_"+GO_str+"_"+cat+"_ensgids.csv",genes_list_tmp,delimiter=", ",fmt="% s")

