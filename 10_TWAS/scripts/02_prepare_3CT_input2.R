# File name: 02_prepare_3CT_input2.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-28 

library(data.table)
library(tidyverse)

MYWD="./10_TWAS/data"
setwd(MYWD)

#- GWAS
df.gwas0 <- fread("gwas_sumstats_scz2022.txt") 
colnames(df.gwas0) <- c("SNP","chr","BP","A1","A2","beta","se")

#- HRC
df.hrc0 <- fread("./workflow/HRC.chr.bp.rsID.gz")

#- input 2
for (i in c(1:22)){
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
  
  a0 <- paste(MYWD,"3CT/eQTL_p0.1_SNP_union/3CT_union_chr",i,".eQTL_combined_p0.1.txt", sep="")
  df.union.snp <- fread(a0, header = F) %>%
    left_join(df.hrc, by=c("V1"="key")) %>% # map for SNP rsID
    filter(is.na(SNP)==F) # remove those without rsID
  # restrict GWAS SNPs to those in the union set of eQTL p<0.1 SNPs
  df.gwas.final <- df.gwas %>%
    filter(SNP %in% df.union.snp$SNP)
  # output per chr gwas results
  outpath2 <- paste(MYWD,"3CT/GWAS/INPUT2_chr",i,"_gwas_scz2022.txt", sep="")
  fwrite(df.gwas.final,outpath2,sep="\t")
  # output SNP ID for later use (updating HRC bfiles)
  outpath2a <- paste(MYWD,"snp_list_perCHR_3CT/snp.list_eQTLp0.1_chr",i,sep="")
  fwrite(df.union.snp %>% distinct(SNP), outpath2a, sep="\t", col.names = F)
}