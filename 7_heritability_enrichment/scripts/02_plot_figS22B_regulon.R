# 02_plots_figS22B_regulon.R

library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./7_heritability_enrichment/output/figS22B/")
a <- fread("df_pLDSC_regulon.tsv")

# plot p-value of enrichment, 
p <- ggplot(a, aes(x=reorder(V1,-Enrichment_p), 
                    y=-log10(Enrichment_p))) +
  geom_bar(stat = "identity", aes(fill = if.sig.fdr)) +
  theme_minimal() +
  labs(
    x = "regulons",
    y = expression("SNP-h"^2 * " enrichment, -log"[10]*"(P value)"),
    fill="FDR≤0.05"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# output plot
ggsave(filename="figS22B_p_regulon_pLDSC.pdf",
       plot=p, width = 5, heigh=6)

#-- end --#