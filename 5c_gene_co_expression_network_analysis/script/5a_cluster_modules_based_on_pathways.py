# author: Lisa Bast
# date: 2024-01-12,  13:46:21
# version: 0.0.1
# about: optionally save sign pathways per module and cluster them, plot sankey plot and cluster map
#        plot excitatory L2-3 IT II red module pathways
#        plot bordeaux module pathways
#        Fig. 5 D,E,H S19D

import pandas as pd
import os
import sys
import numpy as np
import seaborn as sns
#import matplotlib.pyplot as plt

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

path_main = path_code.replace("5c_gene_co_expression_network_analysis\\script","")
path_main_project = path_code.replace("script","")

n_cl = 15
opt_save_modules_pathways = True
GO_strings = ['CC','MF','BP']
sel_strings = ["all_Modules"]#,"all_DEGs","overlap_only"]
opt_modules = "scz_implicated"#"all" # 
opt_test = "LR" #"bimod"
opt_order_by_cluster = False
opt_order_by_ct = True

CT_list_modules = ["Oligodendrocytes" ,"Excitatory_Layer_5-6_CT_and_NP_neurons","Inhibitory_SST_neurons","Inhibitory_VIP_neurons","Inhibitory_LAMP5_neurons",
                "Inhibitory_PVALB_neurons","Astrocytes","Oligodendrocyte_progenitor_cells",
                "Excitatory_Layer_2-3_IT_neurons_I","Excitatory_Layer_2-3_IT_neurons_II",       
                "Excitatory_Layer_3-4_IT_neurons","Excitatory_Layer_3-6_IT_neurons","Excitatory_Layer_5-6_IT_neurons_I"]
                #no modules as cell type populations are too small: "Endothelial_and_mural_cells","Excitatory_Layer_5_6_IT_neurons_II","Microglial_cells",
#paths
path_data = path_main + "/4_data_integration_and_cell_type_annotation/output/"
path_modules_raw = path_main_project + "/output/cluster_name_15CTs/scz/cor_bicor/" 

DF_colors = ut.get_CT_colors_as_df(n_cl,path_data)
DF_colors = ut.get_short_CT_names(DF_colors,"celltype")

if opt_save_modules_pathways:
    ut.save_modules_pathways(GO_strings,path_modules_raw,CT_list_modules,opt_test,opt_modules)
    #otherwise directly load

for GO_str in GO_strings:
    #load modules pathways
    DF_modules = pd.read_csv(path_modules_raw+"DF_sign_modules_"+GO_str+"_all_cell_types_"+opt_test+".csv")

    #for some reason many pathways kicked out
    DF_modules = ut.get_short_CT_names(DF_modules,"Celltype_Module")
    DF_modules['Celltype_Module_short'] = DF_modules['Celltype_Module_short'].str.replace('_',' ')
  
    #plot upset for pathway overlap between modules:
    ut.plot_number_sign_pathways_as_upset_modules(DF_modules,opt_test,GO_str,path_modules_raw,1,True)  

    for sel in sel_strings:
        if sel=="all_Modules":
            #read parent term mapping
            df_parent_term_mapping = pd.read_csv(path_modules_raw+"df_parent_term_mapping_"+GO_str+".csv")
            #visualize number pathways per category normalized by number of pathways in this category
            ut.plot_module_pathway_heatmap_and_clustermap(DF_modules, GO_str, df_parent_term_mapping, path_modules_raw, "",opt_test,False,opt_order_by_cluster,opt_order_by_ct)
        for CT_class in ["Excitatory","Inhibitory","Non-Neuronal"]:
            if CT_class == "Non-Neuronal":
                DF_modules_class = DF_modules[np.logical_and(DF_modules["Celltype_Module_short"].str.startswith("Exc")==False,DF_modules["Celltype_Module_short"].str.startswith("Inh")==False)].copy()
            else:
                DF_modules_class = DF_modules[DF_modules["Celltype_Module"].str.startswith(CT_class)].copy()
            #DF_modules_class = fg.get_simplified_parent_terms(DF_modules_class)
            df_modules_only = DF_modules_class.groupby(['Celltype_Module_short',"parentTerm"]).size().reset_index(name='Count').rename(columns={"Celltype_Module_short":"source","parentTerm":"target","Count":"value"})
            DF_colors = ut.get_CT_colors_as_df(n_cl,path_data)
            DF_colors = ut.get_short_CT_names(DF_colors,"celltype")
            ut.plot_sankey_pathways_2_levels(df_modules_only,GO_str,CT_class,DF_colors,True)


## plot excitatory L2-3 IT II red module pathways
#settings:
opt_export_list_of_module_genes = True # False
opt_plot_ExcL23II_GO_terms = True
opt_export_results_ExcL23II_GO_terms = True

#read Go and syngo results
res_file_str = "_genes_Excitatory_Layer_2-3_IT_neurons_II_red"
results_path_ExcL23II = path_modules_raw+"/Excitatory_Layer_2-3_IT_neurons_II/norm_basic/"
results_path_SYNGO = results_path_ExcL23II+"/GO_SYNGO/"

if opt_export_list_of_module_genes:
    #export list of module genes to run SYNGO
    module_genes = pd.read_csv(path_modules_raw+"/module_gene_mapping_info_Excitatory_Layer_2-3_IT_neurons_II.csv")
    module_genes[module_genes["module_name"].str.endswith("red")]["accession"].to_csv(results_path_ExcL23II+"genes_Excitatory_Layer_2-3_IT_neurons_II_red.txt",index=False,header=False)

if opt_plot_ExcL23II_GO_terms:
    go_file = "Excitatory_Layer_2-3_IT_neurons_II/norm_basic/g_profiler_results_Excitatory_Layer_2-3_IT_neurons_II.csv"
    syngo_files = results_path_SYNGO+"/syngo_ontologies_with_annotations_matching_user_input.xlsx"
    GO_term_category_mapping_file = "T9_Term_category_mapping.xlsx"
    #"module_list_rrvgo_red_BP_size_all.csv"
    DF_GO = pd.read_csv(path_modules_raw+go_file) #add category
    DF_GO = DF_GO[DF_GO["source"]=="GO:BP"]
    DF_GO.drop(columns = [c for  c in DF_GO.columns.to_list() if c not in ["term_id","term_name","pval_Excitatory_Layer_2-3_IT_neurons_II_module_red"]],inplace=True)
    DF_mapping = pd.read_excel(path_data+GO_term_category_mapping_file)
    DF = []
    DF = pd.merge(left = DF_GO, right = DF_mapping, left_on="term_name",right_on="term",how="left")
    DF.drop_duplicates(inplace=True)
    #keep only significant:
    DF = DF[DF["pval_Excitatory_Layer_2-3_IT_neurons_II_module_red"]<0.05]
    #plot barplot of p-values or enrichment:
    if len(DF)>0:
        if "negative_log10_of_adjusted_p_value" not in DF.columns:
            DF["negative_log10_of_adjusted_p_value"] = np.log10(DF["pval_Excitatory_Layer_2-3_IT_neurons_II_module_red"])*(-1)
        ut.plot_GOs_as_horizontal_barplot(DF,res_file_str,results_path_ExcL23II,"term_name","negative_log10_of_adjusted_p_value","GO BP",n_top_pws = 30)
    #SYNGO
    DF_SYNGO = pd.read_excel(syngo_files)
    DF_SYNGO = DF_SYNGO[DF_SYNGO["GSEA 'gene cluster' FDR corrected p-value"]<0.05][["GO domain","GO term name","GSEA 'gene cluster' FDR corrected p-value"]]
    DF_S = DF_SYNGO[["GO term name","GSEA 'gene cluster' FDR corrected p-value"]]#.set_index("GO term name")
    DF_S.rename(columns={"GSEA 'gene cluster' FDR corrected p-value":"FDR corrected p-value"},inplace=True)
    DF_S["negative log10 FDR corrected p-value"] = np.log10(DF_S["FDR corrected p-value"])*(-1)
    if len(DF_S)>0:
        ut.plot_GOs_as_horizontal_barplot(DF_S,res_file_str,results_path_ExcL23II,"GO term name","negative log10 FDR corrected p-value","SYNGO", n_top_pws = 999)
        #DF_SYNGO_p = pd.pivot(data=DF_SYNGO,columns="GO term name",values="GSEA 'gene cluster' FDR corrected p-value")
        ut.plot_sign_syngo_pws(DF_S.set_index("GO term name"),results_path_ExcL23II)

if opt_export_results_ExcL23II_GO_terms:
    #save tables
    with pd.ExcelWriter(results_path_ExcL23II+'sign_pathways_'+res_file_str+'.xlsx') as writer:
        DF.to_excel(writer, sheet_name='sign GO_BP pathways with category')
        DF_S.to_excel(writer, sheet_name='sign SYNGO_pathways')
    

## plot bordeaux module go and syngo pathways
#read Go and syngo results
n_overlapping_genes = ["1","2","4"]
opt_export_genes_overlapping = False
if opt_export_genes_overlapping:
    export_variables = ["gene_ids","accession"]
    gene_filter_options = ["no_filtering","intersection_of_all_modules"]

if opt_export_genes_overlapping:
    for ct in CT_list_modules:
        Module_list = ["Astrocytes_module_yellow",
                "Excitatory_Layer_5-6_IT_neurons_I_module_green",
                "Oligodendrocytes_module_red",
                "Inhibitory_LAMP5_neurons_module_red",
                "Inhibitory_VIP_neurons_module_pink",
                "Excitatory_Layer_2-3_IT_neurons_I_module_red",
                "Excitatory_Layer_2-3_IT_neurons_II_module_black", 
                "Inhibitory_PVALB_neurons_module_magenta"] 
    
        DF_module_genes = pd.read_csv(path_modules_raw+"module_gene_mapping_info_"+ct+".csv")
        DF_module_genes_sel = DF_module_genes[DF_module_genes["module_name"].isin(Module_list)]
        if ct==CT_list_modules[0]:
            DF = DF_module_genes_sel.copy()
        else:
            DF = pd.concat([DF,DF_module_genes_sel])
    for export_var in export_variables:
        for opt_filter_genes in gene_filter_options:
            if opt_filter_genes=="intersection_of_all_modules":
                for th in n_overlapping_genes:
                    DF["count"] = 1
                    DF_p = pd.pivot(DF,index=export_var,columns="module_name",values="count" )
                    DF_p = DF_p.fillna(0)
                    #plt.hist(sorted(DF_p.sum(axis=1),reverse=True))
                    genes_overlapping = DF_p[DF_p.sum(axis=1)>=th].reset_index()
                    genes_overlapping[export_var].to_csv(path_main_project+"/data/"+"bordeaux_modules_genes_"+export_var+"_intersection_of_"+str(th)+".csv",index=False,header=False)
                    print("Genes overlapping: "+str(len(genes_overlapping[export_var].tolist())))
            else:
                DF_unique = DF.drop_duplicates(subset=[export_var])
                DF_unique[export_var].to_csv(path_main_project+"/data/"+"bordeaux_modules_genes_"+export_var+".csv",index=False,header=False)
                print("Genes total: "+str(len(DF[export_var].tolist())))


for n_g in n_overlapping_genes:
    go_file = "/bordeaux_modules/GO/gProfiler_hsapiens_intersecting_magenta_module_genes_"+n_g+".csv"
    syngo_files = "/bordeaux_modules/SYNGO/"+n_g+"_genes_overlap/syngo_ontologies_with_annotations_matching_user_input.xlsx"
    GO_term_category_mapping_file = "Term_category_mapping.xlsx"

    DF_GO = pd.read_csv(path_modules_raw+go_file) #add category
    DF_GO = DF_GO[DF_GO["source"]=="GO:BP"]
    DF_mapping = pd.read_excel(path_modules_raw+GO_term_category_mapping_file)
    DF = []
    DF = pd.merge(left = DF_GO, right = DF_mapping, left_on="term_name",right_on="term",how="left")
    DF.drop_duplicates(inplace=True)
    #plot barplot of p-values or enrichment:
    if len(DF)>0:
        ut.plot_GOs_as_horizontal_barplot(DF,"_bordeaux_modules_"+n_g+"_overlapping_genes",path_modules_raw+"/bordeaux_modules/","term_name","negative_log10_of_adjusted_p_value","GO BP")
    #SYNGO
    DF_SYNGO = pd.read_excel(path_modules_raw+syngo_files)
    DF_SYNGO = DF_SYNGO[DF_SYNGO["GSEA 'gene cluster' FDR corrected p-value"]<0.05][["GO domain","GO term name","GSEA 'gene cluster' FDR corrected p-value"]]
    DF_S = DF_SYNGO[["GO term name","GSEA 'gene cluster' FDR corrected p-value"]]#.set_index("GO term name")
    DF_S.rename(columns={"GSEA 'gene cluster' FDR corrected p-value":"FDR corrected p-value"},inplace=True)
    DF_S["negative log10 FDR corrected p-value"] = np.log10(DF_S["FDR corrected p-value"])*(-1)
    if len(DF_S)>0:
        ut.plot_GOs_as_horizontal_barplot(DF_S,"_bordeaux_modules_"+n_g+"_overlapping_genes",path_modules_raw+"/bordeaux_modules/","GO term name","negative log10 FDR corrected p-value","SYNGO")
        #DF_SYNGO_p = pd.pivot(data=DF_SYNGO,columns="GO term name",values="GSEA 'gene cluster' FDR corrected p-value")
        ut.plot_sign_syngo_pws(DF_S.set_index("GO term name"),path_modules_raw+"/bordeaux_modules/")

    #save tables
    with pd.ExcelWriter(path_modules_raw+"/bordeaux_modules/"+'sign_pathways_'+n_g+'_genes_overlapping.xlsx') as writer:
        DF.to_excel(writer, sheet_name='sign GO_BP pathways with category')
        DF_S.to_excel(writer, sheet_name='sign SYNGO_pathways')