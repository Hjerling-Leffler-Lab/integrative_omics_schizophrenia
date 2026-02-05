#!/bin/bash 

#- name: 01a_pLDSC_calculateLD.sh
#- author: Shuyang Yao (shuyang.yao@ki.se)
#- date: 2024-10-08
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

#- step 1. calculate LD
cd ${WD}/$1/
  for f in *.bed
do
echo $f
# intersect SNPs with the top specific gene (bed files) to get SNPs in the annotations.
intersectBed -c -a $path_name$all_snps -b $f > $f".1000genomes.intersect"
awk '{if($5!=0) print $4}' $f".1000genomes.intersect" > $f".1000genomes.intersect.snp"
mkdir $f"_tissue_dir"
rm $f".1000genomes.intersect" # clean environment
cd $f"_tissue_dir"
# format SNP annotation files, using file "fast_match.py" from J. Bryois (https://github.com/jbryois/scRNA_disease, PMID=32341526)
gunzip *
  for j in $path_name$all_annotations/*.annot
do
echo $j
file_name=`basename $j`
perl $path_name/fast_match.py ../$f".1000genomes.intersect.snp" $f $j > $file_name
done
gzip *annot
# calculate LD scores for annotations
for i in {1..22}
do
sbatch -t 1:00:00 -n 1 -o log_$i --wrap="ldsc.py --l2 --bfile $path_name$plink_file/1000G.EUR.QC.$i --ld-wind-cm 1 --print-snps $path_name$hapmap_snps --annot baseline.$i.annot.gz --out baseline.$i"
done
cd ..
rm $f".1000genomes.intersect.snp" # clean environment

