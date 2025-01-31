#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
myct=args[1]
WD=paste("./9_sc_eQTL/output/16CT/", myct, "/", sep="")
setwd(WD)

library(data.table)
library(tidyverse)
# BiocManager::install("qvalue")
library(qvalue)

#- manage path: use current working directory; i.e., first enter the working directory then run the script
# read data
dat <- fread(paste(myct,"_all.chr.eQTL_combined.txt",sep=""),
             stringsAsFactors = F, data.table = F)
colnames(dat) <- c("phe_id",
                   "phe_chr",
                   "phe_from",
                   "phe_to",
                   "phe_strd",
                   "n_var_in_cis",
                   "dist_phe_var",
                   "var_id",
                   "var_chr",
                   "var_from",
                   "var_to",
                   "nom_pval",
                   "r_squared",
                   "slope",
                   "slope_se",
                   "best_hit")
# estimate q-value
qobj <- qvalue(dat$nom_pval)
dat$qvalue <- qobj$qvalue

summary(dat$nom_pval)
summary(dat$qvalue)

# Output data
fwrite(dat, file = paste(myct,"_eQTL_sumstats.tsv",sep=""),
       sep="\t", col.names=T)
fwrite(dat %>% filter(qvalue<=0.05), file = paste(myct,"_all.chr.eQTL_combined_qvalue0.05.tsv",sep=""),
       sep="\t", col.names=T)

#--- script ends ---#