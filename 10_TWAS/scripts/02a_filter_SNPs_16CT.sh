#!/bin/bash

CT=$1
INDIR= # folder with eQTL full results of 16 cell types.
OUTDIR=./10_TWAS/data/16CT/${CT}
rm ${OUTDIR}/*
for i in {1..22}
do
 cat ${INDIR}/${CT}_eQTL_sumstats.tsv | grep "chr$i\s\+" | awk '$12<=0.1{print $8,$1,$14,$15}' OFS='\t' > ${OUTDIR}/${CT}_chr${i}.eQTL_combined_p0.1.txt
done
