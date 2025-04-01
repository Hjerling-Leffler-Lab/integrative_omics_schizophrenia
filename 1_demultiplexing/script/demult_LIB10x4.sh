
#! /bin/bash -l
#SBATCH -A snic2021-22-317
#SBATCH -p node -n 16
#SBATCH -t 24:00:00
#SBATCH -J demultiplexing_LIB10x4

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_data=/1_demultiplexing/data/

path_cr=/tools/cellranger-6.0.1/

${path_cr}cellranger mkfastq --run=${path_data}200325_A00187_0279_BH77YKDSXY/ --csv=samplesID_LIB10x4.csv --output-dir=/1_demultiplexing/output/library_4/


