#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(tidyverse)

CT=args[1]
MYWD="./10_TWAS/data/3CT/"

#- GWAS
a1 <- "./10_TWAS/data/gwas_sumstats_scz2022.txt"
df.gwas0 <- fread(a1) 
colnames(df.gwas0) <- c("SNP","chr","BP","A1","A2","beta","se")

#- HRC
df.hrc0 <- fread("./workflow/HRC.chr.bp.rsID.gz")

for (i in c(1:22)){
  #- eQTL
  a <- paste(MYWD,CT,"/",CT,"_chr",i,".eQTL_combined_p0.1.txt", sep="")
  
  df.eqtl <- fread(a) %>%
    separate(V1, into=c("chr","pos","ref","alt"), sep=":",remove=F) %>%
    mutate(key=paste(chr,pos,sep=":")) 
  
  #- gwas
  df.gwas <- df.gwas0 %>%
    filter(chr==i)
  
  #- hrc rsID and chr:bp data
  df.hrc <- df.hrc0 %>%
    filter(`#CHROM`==i) %>%
    mutate(key=paste("chr",`#CHROM`,":",POS,sep="")) %>%
    rename(SNP=ID) %>%
    distinct(key, SNP) %>%
    filter(is.na(SNP)==F) %>% 
    filter(SNP!=".") %>%
    group_by(key) %>%
    mutate(n=row_number()) %>%
    ungroup()
  table(df.hrc$n)
  df.hrc <- df.hrc %>% filter(n==1) %>% select(-n)
  
  df.merged <- df.eqtl %>%
    left_join(df.hrc, by="key") 
  
  #- take intercept
  tmp1 <- df.gwas %>% distinct(SNP)
  tmp2 <- df.merged %>% distinct(SNP)
  tmp3 <- df.hrc %>% distinct(SNP)
  snp.list <- tmp1 %>% inner_join(tmp2, by="SNP") %>% inner_join(tmp3, by="SNP")
  
  #- file1: eQTL sumstats
  df.eql.final <- df.merged %>%
    select(SNP, V2, V3, V4) %>%
    rename(Gene=V2, beta=V3, se=V4) %>%
    filter(SNP%in%snp.list$SNP)
  #output: eQTL sumstats: SNP, Gene, beta, se (tab delimited)
  outpath1 <- paste(MYWD,CT,"/INPUT1_chr",i,".eQTL_combined_p0.1.txt", sep="")
  fwrite(df.eql.final,outpath1,sep="\t")
  
  #- file3 and file5: filter by chr, and restrict to only overlapping SNPs, do in plink in the next session.
  
  #- file4: SNPs in eQTL sumstats, same format as bim: V1 (chr, only nr), V2 (SNP_ID), V3 (pos), V4 (pos, same as V3), V5 (minor), V6 (major allele)
  df.bim.final <- df.merged %>%
    mutate(chr2=as.integer(gsub("chr","",chr)),
           pos2=pos) %>%
    distinct(chr2, SNP, pos, pos2, alt, ref) %>%
    filter(SNP%in%snp.list$SNP)
  
  outpath4 <- paste(MYWD,CT,"/INPUT4_chr",i,".eQTL_combined_p0.1.bim", sep="")
  fwrite(df.bim.final,outpath4,sep="\t", col.names=F)
  
}
