
#! /bin/bash -l
#SBATCH -A snic2021-22-317
#SBATCH -p node -n 16
#SBATCH -t 40:00:00

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_ref=/1_demultiplexing/tools/

path_cr=/1_demultiplexing/tools/cellranger-6.0.1/

path_fastq=/1_demultiplexing/output/library_4/H77YKDSXY/

iden=$1
${path_cr}cellranger count --id="Counts_${iden}" --transcriptome=${path_ref}refdata-gex-GRCh38-2020-A --fastqs ${path_fastq}outs/fastq_path/ --sample $1 --expect-cells 5000 --include-introns --outputs-dir /2_alignment/output/

