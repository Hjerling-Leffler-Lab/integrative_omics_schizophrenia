#= 01_get_PCA_GRS.sh
#= SY, 2020-12-21
#= Aim: 
#- get genotype PCs using ricopili pipeline: https://sites.google.com/a/broadinstitute.org/ricopili/pca
#- get GRS, using the PT method based in GRSworkflow (https://github.com/luyi0629/GRSworkflow) on a secured server


#================ SCRIPT BEGINS ================#

#= 0. prepare data and GRSworkflow

GRSworkflow=GRScontainer_v6_optimized/ # download the container and put path to the GRSworkflow folder
cd ${GRSworkflow}


#= 1. GRScontainer settings:

#- 1.1. provide sumstats​
cd ${GRSworkflow}/data/sumstats/raw # put GWAS sumstats here
#- GWAS summary statistics
zcat PGC3_SCZ_wave3_public.v2.tsv.gz | awk 'NR>1{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11}' | sed '1i\
CHR\tSNP\tBP\tA1\tA2\tFRQ_U\tINFO\tOR\tSE\tP
'| tr -s ' ' | tr ' ' '\t' > PGC3_SCZ_wave3_public.v2.tsv

# prepare the sumstats_info.v2.txt file for GRSworkflow
vi ./input-definitions/sumstats_info.v2.txt

#- 1.2 genotype data 
GENODIR=YOUR_PATH_WITH_GENOTYPE_BFILES/ # put the post-QC post-imputation bfiles here

# generate genotype PCA
cd ${GENODIR}
# use ricopili pipeline, pcaer function
pcaer --out brn_pca_noref scz_brn_eur_rk-qc1.hg19.ch.fl.bg_rnaseq_sample_autosome 
# output PCs will be in file brn_pca_noref.menv.trans.mds
# put it under ./8_PRS/output/fig1B/

#- 1.3. prepare input definitions​
cd ${GRSworkflow}/input-definitions
vi step2inputs.csv
# paste genotype bfile name there: scz_brn_eur_rk-qc1.hg19.ch.fl.bg_rnaseq_sample_autosome


#- 1.4. Adapt the workflow to the local slurm, in scripts/​
cd ${GRSworkflow}/scripts

vi RunPipeline.sh
#in RunPipeline.sh: Change the parameters of directory path ‘PROJ_DIR’ & ‘INPUT_PATH’, and edit the ‘PHENO_JOBS’ (# of sumstats) & ’SAMPLE_JOBS’ (# of target geno set) if needed​
# Path to the directory for the sumstats_info.v2.txt and sample2inputs.csv files
INPUT_PATH=GRScontainer_v6_optimized/input-definitions # specify your path to the GRSworkflow folder
# Enter the number of phenotypes for step 1 and samples for step 2
PHENO_JOBS=1 # adjust depending on need.
SAMPLE_JOBS=1
PROJ_DIR=GRScontainer_v6_optimized # specify your path to the GRSworkflow folder

#In  step2.sbatch, step3.sbatch: bind-mount geno data [see below]​
-B YOUR_PATH_WITH_GENOTYPE_BFILES/:YOUR_PATH_TO_GRScontainer_v6_optimized/data/geno/raw \

#- 2. Run PRS calculation:
cd ${GRSworkflow}/
bash scripts/RunPipeline.sh


#--- END ---#