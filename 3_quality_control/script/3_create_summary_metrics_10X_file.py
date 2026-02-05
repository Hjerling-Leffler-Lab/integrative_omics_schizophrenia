# -*- coding: utf-8 -*-
#version: 0.2.0
#Created on Fri Oct  9 14:47:39 2020
#@author: Lisa Bast
#extract QC metrics from several sources (10X summary metrics, filtering result, stats based on raw count data) and summarize them in a dataframe/ .csv file (exports table T3: T3_quality_metrics_per_sample.xlsx)


import os
import pandas as pd
import loompy
import glob
import sys
#import pickle5 as pickle

path_script = os.getcwd()
sys.path.append(path_script)
import utils as ut

path_main = path_script.replace("3_quality_control\\script","")

opt_add_neuron_percentage = False#True
library_nr = range(1, 6) #[1]#
n_samples_per_library = [10, 19, 18, 20, 19]#[2] #
samples_to_drop = ['10x2_19']

path_output = path_main+'/3_quality_control/output/'
path_data = path_main+'/2_alignment/output/'
    
opt_concatenate_metrics_summary=False#True
opt_integrate_filtering_results=True
opt_integrate_raw_data_stats=True

#read metrics exported from 10X summary
df = pd.read_csv(path_data+'metrics_summary_concatenated.csv')

# add IDs to df
df_IDs = pd.read_excel(path_data+"T1_basic_donor_information.xlsx", "basic_donor_information", engine='openpyxl') 
df_IDs.rename(columns = {'scz_status':'Group'},inplace=True)
df_IDs.loc[df_IDs["Group"]=='scz']["Group"]="SCZ"
df_IDs.loc[df_IDs["Group"]=='ctrl']["Group"]="CTRL"
df_IDs['Library'] = df_IDs['Internal_donor_batch_ID'].str[:4]
df_metrics = pd.merge(df_IDs[['Internal_donor_batch_ID', 'donor_ID_universal','Internal_donor_ID','Group','Library']],df,on='Internal_donor_batch_ID')
df_metrics['p_unmapped_reads'] = 100-df_metrics['Fraction Reads in Cells'].str[:-1].astype(float)

#rename important columns
df_metrics.rename(columns = {'Estimated Number of Cells':'num_umi', 
                              'Number of Reads':'num_reads', 
                              'Mean Reads per Cell':'mean_reads_per_umi',
                              'Reads Mapped Confidently to Transcriptome': 'p_mapped_reads',
                              'Reads Mapped Confidently to Intergenic Regions': 'p_genome_not_gene',
                              'Valid Barcodes': 'p_valid_barcodes',
                              'Sequencing Saturation': 'p_sequencing_saturation'
                              },inplace=True)

#type conversion, remove % signs etc
df_metrics['p_mapped_reads'] = df_metrics['p_mapped_reads'].str[:-1].astype(float)
df_metrics['p_genome_not_gene'] = df_metrics['p_genome_not_gene'].str[:-1].astype(float)
df_metrics['p_valid_barcodes'] = df_metrics['p_valid_barcodes'].str[:-1].astype(float)
df_metrics['p_sequencing_saturation'] = df_metrics['p_sequencing_saturation'].str[:-1].astype(float)
df_metrics['num_umi'] = df_metrics['num_umi'].str.replace(',', '').astype(int)
df_metrics['mean_reads_per_umi'] = df_metrics['mean_reads_per_umi'].str.replace(',', '').astype(int)
df_metrics['num_reads'] = df_metrics['num_reads'].str.replace(',', '').astype(int)

#delete irrelevant metrics
del df_metrics['Median Genes per Cell']
#del df_metrics['Valid Barcodes']
#del df_metrics['Sequencing Saturation']
del df_metrics['Q30 Bases in Barcode']
del df_metrics['Q30 Bases in RNA Read']
del df_metrics['Q30 Bases in UMI'] 
del df_metrics['Reads Mapped to Genome']
del df_metrics['Reads Mapped Confidently to Genome']
del df_metrics['Reads Mapped Confidently to Intronic Regions']
del df_metrics['Reads Mapped Confidently to Exonic Regions']
del df_metrics['Reads Mapped Antisense to Gene']
del df_metrics['Fraction Reads in Cells']
del df_metrics['Total Genes Detected']
del df_metrics['Median UMI Counts per Cell']

df_metrics.to_csv(path_output+'metrics_summary_tidy.csv',index=False)

#add other metrics from different files
if opt_integrate_filtering_results:
    #% cells removed in filtering
    df_filt_result_pagoda = pd.read_pickle(path_output+"filtered/df_pagoda_filtering_result_cellranger")
    df_filt_result_pagoda.rename(columns = {'Donor ID':'Internal_donor_ID','Percentage Cells discarded':'p_cells_del_filt_pagoda'},inplace=True)
    #to do: calculate number of cells prior to any filtering based on pagoda df
    df_filt_result_pagoda["Number Cells prior filtering"] = ((df_filt_result_pagoda["Number of cells removed (pagoda filtering)"]*100)/df_filt_result_pagoda["p_cells_del_filt_pagoda"])
    #merge with other filtering data frames:
    df_filt_result_TH = pd.read_pickle(path_output+"filtered/df_pagoda_TH_filtering_result_cellranger")
    df_filt_result_TH.rename(columns = {'Donor ID':'Internal_donor_ID','Percentage Cells discarded':'p_cells_del_filt_TH'},inplace=True)
    df_filt_result_D = pd.read_pickle(path_output+"filtered/df_pagoda_TH_and_D_adj_filtering_result_cellranger")
    df_filt_result_D.rename(columns = {'Donor ID':'Internal_donor_ID','Percentage Cells discarded':'p_cells_del_filt_D'},inplace=True)
    df_filt = pd.merge(df_filt_result_pagoda,df_filt_result_TH,on=['Internal_donor_ID',"Disease"],how='left')
    df_filt = pd.merge(df_filt, df_filt_result_D,on=['Internal_donor_ID',"Disease"],how='left')
    # add numbers of cells removed/ kept in each filtering step
    df_filt['Number of cells post pagoda filtering'] = df_filt['Number Cells prior filtering']-df_filt["Number of cells removed (pagoda filtering)"]
    df_filt['Number of cells removed (TH filtering)'] = (df_filt["p_cells_del_filt_TH"]/100)*df_filt["Number of cells post pagoda filtering"]
    df_filt["Number of cells post TH filtering"] = df_filt["Number of cells post pagoda filtering"]-df_filt["Number of cells removed (TH filtering)"]
    df_filt['Number of cells post D filtering'] = df_filt["Number of cells post TH filtering"]-df_filt["Number of cells removed (predicted doublets)"]
    # calculate percentage of cells discarded overall:
    df_filt["p_cells_del_filt"] = (df_filt_result_pagoda["Number Cells prior filtering"]-df_filt['Number of cells post D filtering'])/df_filt_result_pagoda["Number Cells prior filtering"]

    #merge filtering result with other metrics data frame
    df_metrics1 = pd.merge(df_metrics,df_filt[['p_cells_del_filt','Internal_donor_ID']],on='Internal_donor_ID',how='left')

#read patient-level data
df_patientData = pd.read_excel(path_data+'T1_basic_donor_information.xlsx', 'sample_info',engine='openpyxl') 
df_patientData.rename(columns = {'donor_ID_internal_8':'Internal_donor_ID'},inplace=True)
df_metrics = pd.merge(df_metrics1,df_patientData[['gender','age','PMI_h','Internal_donor_ID']],on='Internal_donor_ID',how='left')

#transform PMI
PMI = df_metrics['PMI_h'].astype('str')
PMI_n=PMI
for p_id,p in enumerate(PMI):
    p_split = p.split(':')
    if len(p_split) ==3:
        PMI_n[p_id] = str(int(p_split[0])+int(p_split[1])/60+int(p_split[2])/(60*60))
df_metrics['PMI_h']=PMI_n.astype('float')
        

if opt_integrate_raw_data_stats:
    #std_reads_per_umi
    df_raw_stats = pd.read_pickle(path_output + "raw/df_stats")
    df_raw_stats.rename(columns = {'Mean counts per barcode':'mean_counts_per_barcode','Std counts per barcode':'std_counts_per_barcode','Median genes per cell':'median_gpc','Donor ID':'Internal_donor_ID'},inplace=True)
    df_metrics = pd.merge(df_metrics,df_raw_stats[['mean_counts_per_barcode','std_counts_per_barcode','median_gpc','Internal_donor_ID']],on='Internal_donor_ID',how='left')

if opt_add_neuron_percentage == True:
    os.chdir(path_main + '/4_data_integration_and_cell_type_annotation/output/15_CTs/')
    files_all = glob.glob('Samples_*.loom')
    bool_first = True
    for f in files_all:
        if "_downsampled_" not in f:
            with loompy.connect(f) as D:
                sample_IDs,n_nuclei,_,_ = ut.get_sampleIDS_nCells_diseaseStatus_sexStatus(D,'human')
            if "neurons" not in f:
                n_cells = n_nuclei
                n_neurons = 0*n_nuclei
            else: 
                n_cells = n_nuclei
                n_neurons = n_nuclei
            #create dataframe with sample_ID and n_cells and neuronal status:
            if bool_first:
                df_neuron_percentage = pd.DataFrame(data = {'Internal_donor_ID': sample_IDs, 'n_cells': n_cells, 'n_neurons': n_neurons})
                bool_first=False
            else:
                for s in sample_IDs:
                    df_neuron_percentage['n_cells'][df_neuron_percentage['Internal_donor_ID']==s] = df_neuron_percentage['n_cells'][df_neuron_percentage['Internal_donor_ID']==s].copy()+n_cells[sample_IDs==s][0]
                    df_neuron_percentage['n_neurons'][df_neuron_percentage['Internal_donor_ID']==s] = df_neuron_percentage['n_neurons'][df_neuron_percentage['Internal_donor_ID']==s].copy()+n_neurons[sample_IDs==s][0]
    df_neuron_percentage['p_neurons'] = (df_neuron_percentage['n_neurons']/df_neuron_percentage['n_cells'])*100
    df_metrics_new = pd.merge(df_metrics,df_neuron_percentage[['p_neurons','Internal_donor_ID']],on='Internal_donor_ID',how='left')
    
    #remove row if Internal_donor_ID is nan
    df_metrics_new = df_metrics_new.dropna(subset=['Internal_donor_ID'])
    df_metrics_new.to_excel(path_output+'T3_quality_metrics_per_sample.xlsx',index=False)
else:    
    #remove row if Internal_donor_ID is nan
    df_metrics = df_metrics.dropna(subset=['Internal_donor_ID'])
    df_metrics.to_excel(path_output+'T3_quality_metrics_per_sample.xlsx',index=False)
