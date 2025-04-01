
#! /bin/bash -l
#SBATCH -A snic2021-22-317
#SBATCH -p node -n 16
#SBATCH -t 24:00:00
#SBATCH -J demultiplexing_LIB10x5

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_data=/1_demultiplexing/data/

path_cr=/tools/cellranger-6.0.1/

${path_cr}cellranger mkfastq --run=${path_data}200325_A00621_0201_AH72CHDSXY/ --csv=samplesID_LIB10x5.csv --output-dir=/1_demultiplexing/output/library_5/

