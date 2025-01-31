#!/bin/bash 

#- name: 01b_pLDSC_partition_h2.sh
#- aim: run partitioned LDSC for different annotations.
#- tool: ldsc (https://github.com/bulik/ldsc)
#- specify $1: the pipeline folder for the annotaitons (e.g., Pipeline_DEG, Pipeline_Specificity, Pipeline_modules)
#- prepare: download the necessary files (below) from ldsc webpage to ./workflow/

#Neccessary files
path_name="./workflow/" # path for the following files (downloaded from ldsc), change if needed
all_snps="1000genomes_phase3_SNPs.bed2"
all_annotations="1000G_EUR_Phase3_baseline_annot"
plink_file="1000G_EUR_Phase3_plink"
hapmap_snps="hm_snp.txt"

WD="./7_heritability_enrichment/data" # specify working directory; structure under this folder: GWAS_sumstats, Pipeline_*/
      # under Pipeline_*/: bed files (derived from previous step in the Rmd file), folder Results/

#- step 2. partition heritability
cd ${WD}/$1/
sumstats=$path_name/scz2022.sumstats.gz # processed summstats (processed by ldsc munge)
gwas_name=`basename $sumstats | cut -d "." -f 1`
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
all_annotations="1000G_EUR_Phase3_baseline"

for f in `ls ${WD}/$1/*.bed`
do
echo $f
cd ${WD}/$1/${f}_tissue_dir
sbatch -t 1:00:00 -n 1 --mem=20g -o log_partitioned_h2_$gwas_name --wrap="ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ${WD}/$1/Results/${gwas_name}_${f}"
done
