# File name: rscript_01_b2_Get_Exp_Input_3CT.R
# Aim: generate expression input for eQTL analysis.

#--- script starts ---#

# Libraries and prepare to load loom files
library(tidyverse)
library(data.table)
library("rhdf5")
library("snow")

# specify paths
filepath <- "" # loom files
WD <- "./9_sc_eQTL/data/eQTL/pheno/16CT/"
setwd(WD)
celltype.list <- list.files(WD, pattern = "gene_passExpFilter_15cluster_")
celltype.list <- gsub("gene_passExpFilter_15cluster_","",celltype.list)
celltype.list <- gsub(".tsv","",celltype.list)

# read gene matrix, needed in loop
geneMx0 <- fread("./workflow/geneMatrix.tsv",
                 stringsAsFactors=F,
                 data.table=F) %>%
  filter(gene_type=="protein_coding",
         is.na(hg19g0)==FALSE,
         is.na(g1)==FALSE,
         is.na(g2)==FALSE,
         PAR=="FALSE",
         hg19g0!="chrM") %>%
  mutate(tss=as.integer(ifelse(gstr=="+", g1, g2)), # 0-based
         tss0=as.integer(tss-1)) %>%
  select(ensgid, hg19g0, tss0, tss, gstr)

# read sample info, needed in loop
sample.df <- fread("./workflow/brn_GRS_EUR_donorIDinternal8.tsv", # get through .R script ./9_sc_eQTL/scripts/99_brn_GRS_EUR_donorIDinternal8.R
                   stringsAsFactors = F, data.table=F) 

#- start loop
for (i in 1:length(celltype.list)){
  myCT = celltype.list[i]
  # 1. gene filters
  # read geneMatrix, filter for pr-coding
  # read cell type gene lists that passed expression filter
  # intersect
  genelist <- fread(paste(WD,"gene_passExpFilter_15cluster_",myCT,".tsv", sep=""),
                    stringsAsFactors=F, data.table=F)
  # 11932 genes
  geneMx <- geneMx0 %>% filter(ensgid %in% genelist$ENSGID)
  nrow(geneMx)
  
  # 2. Get expression
  # loop over loom files of the same type
  # per loom file: 
  #- filter genes
  #- aggregate per donor
  #- make to long format
  #- merge with the previous file and aggregate per donor
  #- continue loop
  
  # myCT
  filelist <- grep(list.files(filepath, pattern = paste(myCT,".loom",sep="")), pattern='downsampled', invert=TRUE, value=TRUE)
  file=paste(filepath, filelist, sep="")
  h5f <- H5Fopen(file, flags="H5F_ACC_RDONLY")
  cell_donor <- cbind(h5f$col_attrs$CellID,
                      h5f$col_attrs$Donor) %>%
    as.data.frame()
  colnames(cell_donor) <- c("CellID", "Donor")
  
  exp <- as.data.frame(t(h5f$matrix))
  colnames(exp) <- h5f$col_attrs$CellID
  exp$Accession <- h5f$row_attrs$Accession
  exp <- exp %>% 
    filter(Accession %in% geneMx$ensgid) %>%
    select(Accession, everything())
  exp.long <- exp %>% 
    pivot_longer(col=-Accession, 
                 names_to = "CellID",
                 values_to = "count")
  rm(exp)
  exp.long <- exp.long %>% 
    left_join(cell_donor, by="CellID") %>% 
    group_by(Accession, Donor) %>%
    summarise(agg_count=sum(count)) %>%
    ungroup()
  colnames(exp.long)[1] <- "ENSGID"
  # output aggregated data
  outfile <- paste(WD,myCT,"_agg_gene_donor_count.tsv", sep="")
  fwrite(exp.long, outfile, sep="\t", col.name=T)
  
  # 3. Format
  # convert to tmp per donor 
  exp.long <- exp.long %>%
    group_by(Donor) %>%
    mutate(total_count=sum(agg_count))%>%
    mutate(exp_tpm=agg_count*1000000/total_count) %>%
    ungroup()
  exp.long <- exp.long %>%
    left_join(sample.df, by=c("Donor"="donor_ID_internal_8"))
  
  # fix sample order to be the same with the order in covariate and genotype files
  exp.long <- exp.long %>%
    filter(is.na(IID)==F) %>%
    mutate(batch=ifelse(substr(IID,4,4)=="-","batch1", "batch2")) %>%
    arrange(batch, IID, ENSGID)
  exp.long %>% distinct(IID) %>% nrow() # n=73
  
  df <- exp.long %>%
    filter(is.na(IID)==F,
           IID!="") %>% # filter out rows: those without IID (e.g., non-EUR)
    select(ENSGID, IID, exp_tpm) %>% # filter for useful columns
    pivot_wider(names_from = IID,
                values_from = exp_tpm) # to wide format, will be columns 4 and above in 
  df <-  df %>%
    inner_join(geneMx, by=c("ENSGID"="ensgid")) %>%
    mutate(ENSGID1=ENSGID) %>%
    select(hg19g0, tss0, tss, ENSGID1, ENSGID, gstr, everything()) %>%
    arrange(hg19g0, tss0, tss)
  colnames(df)[1:6] <- c("#Chr", "start", "end", "pid", "gid", "strand")
  
  # save per chr for eQTL runs.
  clusteroutpath=paste(WD,myCT,sep="")
  system(paste("mkdir ",clusteroutpath,sep=""))
  for (j in 1:22){
    chrfilename=paste(clusteroutpath,"/",myCT,"_exp_tpm.chr",j,"_TSS.bed",sep="")
    df %>% filter(`#Chr`==paste0("chr",j)) %>%
      fwrite(chrfilename, sep="\t",col.names = T)
  }
}

#--- end ---#