#!/bin/bash

i=$1 

# get cell type list
CelltypeList=(Astrocytes Endothelial_and_mural_cells Excitatory_Layer_2-3_IT_neurons_I Excitatory_Layer_2-3_IT_neurons_II Excitatory_Layer_3-4_IT_neurons Excitatory_Layer_3-6_IT_neurons Excitatory_Layer_5-6_CT_and_NP_neurons Excitatory_Layer_5-6_IT_neurons_I Excitatory_Layer_5-6_IT_neurons_II Inhibitory_LAMP5_neurons Inhibitory_PVALB_neurons Inhibitory_SST_neurons Inhibitory_VIP_neurons Microglial_cells Oligodendrocyte_progenitor_cells Oligodendrocytes)

# working dir
cd ./10_TWAS/data/16CT/

# get the union of SNPs per CT per chr
for CT in $CelltypeList
do
 echo ${CT}
 mkdir ${CT}
 cat ${CT}/${CT}_chr${i}.eQTL_combined_p0.1.txt | awk '{print $1}' >> eQTL_p0.1_SNP_union/tmp_16CT_union_chr${i}.eQTL_combined_p0.1.txt 
done