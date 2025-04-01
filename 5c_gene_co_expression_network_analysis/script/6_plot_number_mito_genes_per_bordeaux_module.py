# author: Lisa Bast
# date: 2024-11-04,  17:11:39
# version: 0.0.1
# about: plot number of mitochondrial DNA genes and number of nuclear mitochondrial genes per bordeaux module
#        export genes
#        Fig. S22A

import os
import sys

## specify project path:
os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
sys.path.append(path_code)
import utils as ut

#paths:
path_project_main = path_code.replace("script","")
path_main = path_code.replace("6_cross_and_inter_modality_comparisons\\script","")

path_GO_results_BM = path_project_main + '/output/cluster_name_15CTs/scz/cor_bicor/variable/bordeaux_modules/'
path_module_genes = path_main + "/5c_gene_co_expression_network_analysis/output/cluster_name_15CTs/scz/cor_bicor/"
path_genes_matrix = path_project_main + "/output/Gene_lists/GO_DB/"

DF_BM_genes = ut.get_bordeaux_module_genes_with_pathways_per_category(path_GO_results_BM,path_module_genes,path_genes_matrix)
#filter for category Mitchondria/ ATP synthesis:
DF_BM_genes_mito = DF_BM_genes[DF_BM_genes["category"]=="Mitochondria/ ATP synthesis"]
#get short module names
DF_BM_genes_mito = ut.get_short_CT_names(DF_BM_genes_mito,"Celltype_Module")
ut.plot_number_nuclear_mito_and_mito_genes_per_bordeaux_module(DF_BM_genes_mito[["gene_name","Celltype_Module_short"]],path_project_main+"/output/")
all_mito_pw_genes = DF_BM_genes_mito["gene_name"].unique()

all_mito_pw_genes.tofile(path_project_main+"/output/"+"all_nuclear_and_mitochondrial_genes_of_bordeaux_modules.csv",sep=",")