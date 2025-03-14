# -*- coding: utf-8 -*-
# extract 10X summary metrics from each library and sample and save in one file
# date: Fri Oct  9 12:48:39 2020
# @author: Lisa Bast
# version: 0.3.2

import os
import pandas as pd

path_project = path_code.replace("2_alignment\\script","")

library_nr = range(1, 6) #[1]#
n_samples_per_library = [10, 19, 18, 20, 19]#[2] #
samples_to_drop = ['10x2_19']
path1 = path_project + '/2_alignment/output/'


bool_first_file = True
for l_id, lib in enumerate(library_nr):
    for s_id in range(1,n_samples_per_library[l_id]+1):
        os.chdir(path1+'Library_'+str(lib)+'/Counts_10x'+str(lib)+'_'+str(s_id)+'/outs/')
        df = pd.read_csv('metrics_summary.csv')
        df['donor_ID_library']='10x'+str(lib)+'_'+str(s_id)
        if bool_first_file == True:
            df_all = df.copy()
            bool_first_file = False
        else:
            df_all = pd.concat([df_all, df])

os.chdir(path1)
df_all.to_csv(path1+'metrics_summary_concatenated.csv',index=False)
 