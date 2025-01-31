# author: Lisa Bast
# date: 2024-05-28,  11:43:02
# version: 0.0.2
# about: plot age distribution, Fig. 1A

import pandas as pd
import seaborn as sns
import os
import sys

def switch_time_format(str_time):
    if ":" in str_time: 
        str_array = str_time.split(':')
        new_str = int(str_array[0])+int(str_array[1])/60+int(str_array[2])/(60*60)
    else:
        new_str = str_time
    return new_str

os.chdir(os.path.dirname(__file__))
path_code = os.getcwd()
os.chdir("../")
path_preprocessing = os.getcwd()
os.chdir("../")
path_main = os.getcwd()
sys.path.append(path_code)
import utils as ut

patient_data = path_main+"/2_alignment/output/T1_basic_donor_information.xlsx"
path_results = path_preprocessing+"/output/"
df = pd.read_excel(patient_data,sheet_name="basic_donor_information")
df_sel = df[df["scRNAseq_performed"]=="yes"][["age","gender","scz_status","PMI_h"]]

#convert PMI_h to float
df_sel["PMI_h"] = [switch_time_format(st) for st in df_sel["PMI_h"].tolist()]
df_sel["PMI_h"] = df_sel["PMI_h"].astype("float")

ut.plot_distributions(df_sel,path_results,"age")
ut.plot_distributions(df_sel,path_results,"PMI_h")
ut.plot_distributions(df_sel,path_results,"gender")