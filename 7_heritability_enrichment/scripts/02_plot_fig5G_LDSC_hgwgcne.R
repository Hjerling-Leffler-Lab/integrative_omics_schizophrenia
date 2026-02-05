# 02_plot_fig5G_LDSC_hgwgcne.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-06-16
# aim: plot figure 5G

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

ggsave(p1,
       filename="Fig5G.LDSC_enrichment.pdf",
       width = 4, height = 4)

#--- end ---#