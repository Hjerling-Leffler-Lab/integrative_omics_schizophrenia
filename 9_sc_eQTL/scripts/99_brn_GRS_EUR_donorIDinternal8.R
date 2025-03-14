#- brn_GRS_EUR_donorIDinternal8.R
#- shuyang yao, 20211222
#- aim: get IID and donor_ID_internal_8 (used in expression data processing)
#- for only EUR sample, which will be used for eQTL analysis
#- output: brn_GRS_EUR_donorIDinternal8.tsv

library(data.table)
library(tidyverse)
library(readxl)
# read EUR-only IIDs for both batches
df <- fread("./8a_PRS//brn_GRS_EUR.tsv",
            stringsAsFactor=F,
            data.table=F)
df1 <- df %>% filter(batch=="batch1")
df2 <- df %>% filter(batch=="batch2")
# read sample info master excel sheet
sample.xl <- read_excel("./workflow/T1_basic_donor_information.xlsx") %>%
  select(donor_ID_universal, donor_ID_internal_6, donor_ID_internal_8)

# if batch1, merge with donor_ID_internal_6
# if batch2, merge with donor_ID_universal
df1 <- df1 %>% 
  left_join(sample.xl, by=c("IID"="donor_ID_internal_6")) %>%
  select(IID, donor_ID_internal_8)

df2 <- df2 %>% 
  left_join(sample.xl, by=c("IID"="donor_ID_universal")) %>%
  select(IID, donor_ID_internal_8)

df <- rbind(df1, df2) %>%
  arrange(donor_ID_internal_8)
fwrite(df, 
       "./workflow/brn_GRS_EUR_donorIDinternal8.tsv",
       col.names = T, sep="\t")
