# 04a_LDSC_2ModuleBlocks_enrichment.R
library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./7_heritability_enrichment/output/fig5G/")
df <- fread("df_03_pLDSC_29modules_with2moduleBlocks_shortnames.tsv")

p1 <- ggplot(df %>% filter(`Prop._SNPs`>0.01),
              aes(x=reorder(module_shortnames,-order), y=Enrichment, color=as.character(if.EnrichmentFDR.05)))+
  geom_point() +
  geom_errorbar(aes(ymax=(Enrichment+1.96*Enrichment_std_error), 
                    ymin=(Enrichment-1.96*Enrichment_std_error)), width=.2) +
  coord_flip() +
  scale_color_manual(values=c("grey","salmon"), name="FDR sig.", labels=c("no","yes")) +
  geom_hline(yintercept = 1, linetype="dashed", color="grey") +
  theme_classic() +
  xlab("") + ylab("LDSC enrichment (95% CI)")

p2 <- ggplot(df,
              aes(x=reorder(module_shortnames,-order), y=, `Prop._SNPs`, fill=as.character(if.EnrichmentFDR.05)))+
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("grey","salmon"), name="FDR sig.", labels=c("no","yes")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 0.01, linetype="dashed", color="black") 

ggsave(p1,
       filename="Fig5e.LDSC_enrichment.pdf",
       width = 4, height = 4)
ggsave(p2,
       filename="FigS20a.LDSC_module_size.pdf",
       width = 4, height = 4)
