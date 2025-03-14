# author: Lisa Bast
# date: 2024-06-12,  15:32:15
# version: 0.0.1
# about: plot summary graphic for pathway analyses
# to do: add Ruzicka and Batiuk results as columns

import sys 
import pandas as pd
import numpy as np
import math
import os

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()

#load helper functions:
sys.path.append(path_code)
import utils as ut

#specify remaining paths:
path_project_main = path_code.replace("script","")
path_main = path_code.replace("6_cross_and_inter_modality_comparisons\\script","")
path_DEG_result = path_main + "/5b_differential_gene_expression_analysis/output/DEGs/16_CTs/design_with_6_RUVs/"
DEG_filename = "df_DESeq2_all_genes_all_CTs.csv"
term_category_mapping_path = path_main + "/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/"
results_path = path_project_main + "/output/"
data_path_Ruzicka = path_project_main + "/data/Gene_lists/Ruzicka_et_al_2024/"
data_path_Batiuk = path_project_main + "/data/Gene_lists/Batiuk_et_al_2022/"

#settings:
GO_strings = ["BP","CC","MF"]
opt_load_categories = True
opt_add_other_snRNAseq_studies = True

#load data frames:
DF_pw_result = pd.read_csv(results_path+"pathway_results_all_analyses.csv")
for GO_str in GO_strings:
    filename = "reduced_terms_rrvgo_all_"+GO_str+"_size.csv"
    if opt_load_categories and GO_str == "BP":
        DF_rrvgo_cat = pd.read_excel(results_path+"reduced_terms_rrvgo_all_BP_size_with_categories_LB_curated.xlsx")
        DF_rrvgo_cat["category_incomplete"] = DF_rrvgo_cat["category"]
        DF_rrvgo_cat.loc[DF_rrvgo_cat["category"].isna(),"category"] = DF_rrvgo_cat.loc[DF_rrvgo_cat["category"].isna(),"suggestion_Lisa"]
    else:
        DF_rrvgo_cat = pd.read_csv(results_path+filename)
        if GO_str == "BP":
            DF_term_cat_mapping = pd.read_excel(term_category_mapping_path + "Term_category_mapping.xlsx")
            DF_rrvgo_cat = pd.merge(left=DF_rrvgo_cat,right=DF_term_cat_mapping,on="term",how="left")
            DF_rrvgo_cat.to_csv(results_path+"reduced_terms_rrvgo_all_"+GO_str+"_size_with_categories.csv")
    #filter DF_pw_result for respective GO category
    DF_pw = DF_pw_result[DF_pw_result["subgroup"]=="GO:"+GO_str]
    DF_pw_w_cat = pd.merge(left=DF_pw,right=DF_rrvgo_cat,left_on="ID",right_on="go", how="left")
    #count occurences per category:
    DF_pw_w_cat["Number pathways"] = 1
    DF_pw_w_cat.loc[DF_pw_w_cat["TestVar"]=="modules","TestVar"]="any module"
    DF_pw_w_cat.loc[DF_pw_w_cat["TestVar"]=="bordeaux_modules","TestVar"]="any bordeaux module"
    if GO_str=="BP":
        index_var="category"
    else:
        index_var="parentTerm"
    #add other snRNAseq studies:
    if GO_str=="BP" and opt_add_other_snRNAseq_studies:
        #integrate Ruzicka and Batiuk pathways
        pws_Batiuk = pd.read_excel(data_path_Batiuk + "Supplementary_Dataset_Tables_1-4.xlsx",sheet_name = "Table 2",skiprows=2) 
        pws_Batiuk_upregulated = pws_Batiuk[pws_Batiuk["qvalue"]<0.05][["ID","Description","qvalue"]] #higher CT resolution
        pws_Batiuk_downregulated = pws_Batiuk[pws_Batiuk["qvalue.1"]<0.05][["ID.1","Description.1","qvalue.1"]] #higher CT resolution
        pws_Batiuk_upregulated["TestVar"]="DEGs_up_Batiuk"
        pws_Batiuk_downregulated["TestVar"]="DEGs_down_Batiuk"
        pws_Batiuk_upregulated.rename(columns={"Description":"term","qvalue":"P.bonf.group"},inplace=True)
        pws_Batiuk_downregulated.rename(columns={"ID.1":"ID","Description.1":"term","qvalue.1":"P.bonf.group"},inplace=True)
        #merge with other pathway results:
        DF_pw_w_cat = pd.concat([DF_pw_w_cat,pws_Batiuk_upregulated])
        DF_pw_w_cat = pd.concat([DF_pw_w_cat,pws_Batiuk_downregulated])

        pws_Ruzicka = pd.read_excel(data_path_Ruzicka + "science.adg5136_data_s7.xlsx",sheet_name = "Combined_Pathways_Annotated_01-")
        pws_Ruzicka = pws_Ruzicka[pws_Ruzicka["source"]=="GO:"+GO_str]
        pws_Ruzicka_up = pws_Ruzicka[pws_Ruzicka["Direction"]=="Up"][["term_name","ClusterName","p_values"]]
        pws_Ruzicka_down = pws_Ruzicka[pws_Ruzicka["Direction"]=="Down"][["term_name","ClusterName","p_values"]]
        pws_Ruzicka_up.rename(columns={"term_name":"term","ClusterName":"Ruzicka_category","p_values":"P.bonf.group"},inplace=True)
        pws_Ruzicka_down.rename(columns={"term_name":"term","ClusterName":"Ruzicka_category","p_values":"P.bonf.group"},inplace=True)
        pws_Ruzicka_up["TestVar"] = "DEGs_up_Ruzicka"
        pws_Ruzicka_down["TestVar"] = "DEGs_down_Ruzicka"
        pws_Ruzicka_up["subgroup"] ="GO_BP"
        pws_Ruzicka_down["subgroup"] ="GO_BP"
        #merge with other pathway results
        DF_pw_w_cat = pd.concat([DF_pw_w_cat,pws_Ruzicka_up])
        DF_pw_w_cat = pd.concat([DF_pw_w_cat,pws_Ruzicka_down])
        #to do: add GO ids or manually curate category name
        #search for term in same df and if avail, copy ID and category and suggestion_Lisa

        #make term complete:
        for row in DF_pw_w_cat.index[DF_pw_w_cat["term"].isna()].tolist():
            DF_pw_w_cat.loc[row,"term"] = DF_pw_w_cat.iloc[row]["geneset"].split(" ",1)[1].replace('_',' ')
        #copy GO ID if missing
        for row in DF_pw_w_cat.index[DF_pw_w_cat["go"].isna()].tolist():
            #go missing:
            if pd.isnull(DF_pw_w_cat.iloc[row]["ID"])==False:
                DF_pw_w_cat.loc[row,"go"] = DF_pw_w_cat.iloc[row]["ID"]
            else:
                current_term = DF_pw_w_cat.iloc[row]["term"]
                go = DF_pw_w_cat[DF_pw_w_cat["term"]==current_term]["go"].unique().tolist()
                if len(go)==1:
                    DF_pw_w_cat.loc[row,"go"] = go[0]
                elif len(go)>1:
                    print(go)
                    print(current_term)
                    print("------------")
            #ID missing, but not go:
            if pd.isnull(DF_pw_w_cat.iloc[row]["go"])==False and pd.isnull(DF_pw_w_cat.iloc[row]["ID"])==True:
                DF_pw_w_cat.loc[row,"ID"] = DF_pw_w_cat.iloc[row]["go"]
            else:
                current_term = DF_pw_w_cat.iloc[row]["term"]
                ID = DF_pw_w_cat[DF_pw_w_cat["term"]==current_term]["ID"].unique().tolist()
                if len(ID)==1:
                    DF_pw_w_cat.loc[row,"ID"] = ID[0]
                elif len(ID)>1:
                    print(ID)
                    print(current_term)
                    print("------------")

        #copy category from same term
        for row in DF_pw_w_cat.index[np.logical_or(DF_pw_w_cat["suggestion_Lisa"].isna(),DF_pw_w_cat["category"].isna())].tolist():
            current_ID = DF_pw_w_cat.iloc[row]["ID"]
            #find this term elsewhere in the dataframe
            DF_sel = DF_pw_w_cat[np.logical_and(DF_pw_w_cat["ID"]==current_ID,DF_pw_w_cat["category"].isna()==False)]
            category = DF_sel["category"].unique().tolist()
            if len(category)==1:
                DF_pw_w_cat.loc[row,"category"] = category[0]
            elif len(category)>1:
                print(category)
                print(current_term)
                print("------------")

        
        if opt_add_other_snRNAseq_studies:
            #export DF
            DF_pw_w_cat.to_csv(results_path+"reduced_terms_rrvgo_all_"+GO_str+"_size_with_categories_with_other_snRNAseq_studies.csv")

    #make DF pivot
    DF_pw_pivot = pd.pivot_table(data = DF_pw_w_cat,index=index_var,columns = "TestVar",values="Number pathways",aggfunc = "sum")
    #plot heatmaps:
    ut.plot_integrated_analysis_maps(DF_pw_pivot[["longread","TWAS","DEGs_up","DEGs_down","any module","any bordeaux module","proteomics"]],GO_str,results_path)

        