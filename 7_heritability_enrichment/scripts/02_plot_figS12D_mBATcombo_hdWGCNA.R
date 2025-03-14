library(data.table)
library(ggplot2)

setwd("./7_heritability_enrichment/output/figS12D/")

df <- fread("df_01_phyper_mBATenrich_DEG0.3_GO.tsv")
p1 <- ggplot(df, 
             aes(x=reorder(GO_term,-fdr_p_hg), 
                 y=-log10(p_hg), alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  xlab("") + ylab(expression('-log'[10]*'(P)'))
p2 <- ggplot(df, 
             aes(x=reorder(GO_term,-fdr_p_hg), 
                 y=k_nr_genes_sample, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  xlab("") + ylab("Nr. genes")

# output
ggsave(p1, "figS12D_mBATcombo_sig.pdf", width=7, height=7)
ggsave(p2, "figS12D_mBATcombo_nGene.pdf", width=7, height=7)
