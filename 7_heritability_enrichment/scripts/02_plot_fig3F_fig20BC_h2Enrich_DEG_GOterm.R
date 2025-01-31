library(data.table)
library(tidyverse)
library(rstatix)
library(ggplot2)

setwd("./7_heritability_enrichment/output/fig3F/") 

a <- fread("01_annot_size.tsv")
b <- fread("02_LDSC_res_DEG_GO_terms_size0.01.tsv")

p1 <- ggplot(a, aes(x=reorder(annot,`Prop._SNPs`), y=Prop._SNPs)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  xlab("GO terms") +
  geom_hline(yintercept = 0.01, linetype="dashed", color="grey")

p2 <- ggplot(b, aes(x=reorder(annot,-Enrichment_p), y=-log10(Enrichment_p), alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_classic() +
  xlab("GO terms")

p3 <- ggplot(b, aes(x=reorder(annot,-Enrichment_p), y=Enrichment, alpha=if.sig.fdr)) +
  geom_point() +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, 
                    ymax=Enrichment+Enrichment_std_error*1.96), 
                width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  theme_classic() +
  xlab("GO terms") + ylab("Enrichment (95% CI)") +
  geom_hline(yintercept = 1, linetype="dashed", color="grey")

ggsave(filename=paste(outdir,"figS12B_annot_size.pdf"),
       plot=p1, 
       width=7, height=5, dpi = 300)
ggsave(filename=paste(outdir,"fisS12C_Enrichment_p.pdf"),
       plot=p2, 
       width=7, height=2, dpi = 300)
ggsave(filename=paste(outdir,"fig3F_Enrichment_CI.pdf"),
       plot=p3, 
       width=7, height=2, dpi = 300)

