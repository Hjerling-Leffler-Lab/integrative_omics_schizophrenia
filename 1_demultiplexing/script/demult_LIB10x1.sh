#! /bin/bash -l
#SBATCH -A snic2021-22-317
#SBATCH -p node -n 16
#SBATCH -t 24:00:00
#SBATCH -J demultiplexing_LIB10x1

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_data=/1_demultiplexing/data/

path_cr=/tools/cellranger-6.0.1/

${path_cr}cellranger mkfastq --run=${path_data}190712_A00187_0179_AHLJTVDSXX/ --csv=samplesID_LIB10x1.csv --output-dir=/1_demultiplexing/output/library_1/
