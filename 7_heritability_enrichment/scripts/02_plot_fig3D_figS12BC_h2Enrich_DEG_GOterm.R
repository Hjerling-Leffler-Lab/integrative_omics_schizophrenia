# 02_plot_fig3D_figS12BC_h2Enrich_DEG_GOterm.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2025-06-28
# aim: plot Fig 3D and Fig S12B and S12C

library(data.table)
library(tidyverse)
library(rstatix)
library(ggplot2)

setwd("./7_heritability_enrichment/output/fig3D/") 

a <- fread("02_pLDSC_results_DEG_GO_gt05pctSNP_alpha1.tsv")

a2 <- a %>%
  adjust_pvalue(p.col="Enrichment_p", output.col = "Enrichment_fdr_onlyBP", method="fdr") %>%
  mutate(if.sig.fdr.BP=ifelse(Enrichment_fdr_onlyBP<=0.05,"yes","no"))

a3 <- a %>% filter(substr(GO,1,2)=="BP") %>%
  adjust_pvalue(p.col="Enrichment_p", output.col = "Enrichment_fdr_onlyBP", method="fdr") %>%
  mutate(if.sig.fdr.BP=ifelse(Enrichment_fdr_onlyBP<=0.05,"yes","no"))
a3$GO <- gsub("BP ","",a3$GO)

p1 <- ggplot(a2, aes(x=reorder(GO,`Prop._SNPs`), y=Prop._SNPs)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  xlab("GO terms") +
  geom_hline(yintercept = 0.005, linetype="dashed", color="grey")

p2 <- ggplot(a3, aes(x=reorder(GO,-Enrichment_p), y=-log10(Enrichment_p), alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  xlab("GO-BP terms") + ylab("-log10(P),\nSignificance of\nSchizophrenia SNP-h2")

p3 <- ggplot(a2, aes(x=reorder(GO,-Enrichment_p), y=-log10(Enrichment_p), alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  xlab("GO terms") + ylab("-log10(p-value)")

ggsave(filename=paste(WD,"figS12B_annot_size.pdf", sep=""),
       plot=p1, 
       width=4, height=4, dpi = 300)
ggsave(filename=paste(WD,"fig3D_significance_gobp.pdf", sep=""),
       plot=p2, 
       width=5, height=3, dpi = 300)
ggsave(filename=paste(WD,"figS12C_significance_all.pdf", sep=""),
       plot=p3, 
       width=5, height=4, dpi = 300)

#--- end ---#