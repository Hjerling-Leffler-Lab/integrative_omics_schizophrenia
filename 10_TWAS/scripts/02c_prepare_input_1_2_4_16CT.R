library(data.table)
library(tidyverse)

CT=commandArgs(trailingOnly = TRUE)
MYWD="./10_TWAS/data/"

#- GWAS
df.gwas0 <- fread(gwas_sumstats_scz2022.txt) 
colnames(df.gwas0) <- c("SNP","chr","BP","A1","A2","beta","se")

#- HRC
df.hrc0 <- fread("./workflow/REF/HRC.chr.bp.rsID.gz") # folder with HRC data

# file 2: SCZ GWAS--run once for 16CT
if(CT=="Excitatory_Layer_2-3_IT_neurons_I"){
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
    
    a0 <- paste(MYWD,"16CT/eQTL_p0.1_SNP_union/16CT_union_chr",i,".eQTL_combined_p0.1.txt", sep="")
    df.union.snp <- fread(a0, header = F) %>%
      left_join(df.hrc, by=c("V1"="key")) %>% # map for SNP rsID
      filter(is.na(SNP)==F) # remove those without rsID
    # restrict GWAS SNPs to those in the union set of eQTL p<0.1 SNPs
    df.gwas.final <- df.gwas %>%
      filter(SNP %in% df.union.snp$SNP)
    # output per chr gwas results
    #system(paste("mkdir ", MYWD, "16CT/GWAS/",sep=""))
    outpath2 <- paste(MYWD,"16CT/GWAS/INPUT2_chr",i,"_gwas_scz2022.txt", sep="")
    fwrite(df.gwas.final,outpath2,sep="\t")
    # output SNP ID for later use (restricing HRC bfiles)
    outpath2a <- paste(MYWD,"/snp_list_perCHR_16CT/snp.list_eQTLp0.1_chr",i,sep="")
    fwrite(df.union.snp %>% distinct(SNP), outpath2a, sep="\t", col.names = F)
  }
}else{print("GWAS file processed already.")}

#- input files 1 and 4 (eQTL): get SNP rsID for the filtered (p<0.01) SNPs.
for (i in c(1:22)){
  #- eQTL
  a <- paste(MYWD,"16CT/",CT,"/",CT,"_chr",i,".eQTL_combined_p0.1.txt", sep="")
  
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
  outpath1 <- paste(MYWD,"16CT/",CT,"/INPUT1_chr",i,".eQTL_combined_p0.1.txt", sep="")
  fwrite(df.eql.final,outpath1,sep="\t")
  print("File 1 done.")
  
  #- file4: SNPs in eQTL sumstats, same format as bim: V1 (chr, only nr), V2 (SNP_ID), V3 (pos), V4 (pos, same as V3), V5 (minor), V6 (major allele)
  df.bim.final <- df.merged %>%
    mutate(chr2=as.integer(gsub("chr","",chr)),
           pos2=pos) %>%
    distinct(chr2, SNP, pos, pos2, alt, ref) %>%
    filter(SNP%in%snp.list$SNP)
  
  outpath4 <- paste(MYWD,"16CT/",CT,"/INPUT4_chr",i,".eQTL_combined_p0.1.bim", sep="")
  fwrite(df.bim.final,outpath4,sep="\t", col.names=F)
  print("File 4 done.")
}

#-- end --#
