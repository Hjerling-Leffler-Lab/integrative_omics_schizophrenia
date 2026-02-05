# 02_plot_fig2E_h2SpecificGenes.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-07
# aim: plot figure 2E

library(data.table)
library(tidyverse)
library(ggplot2)

# path
mypath <- "./7_heritability_enrichment/output/fig2E/"

e <- fread(paste(mypath,"df_h2SpecificGenes.tsv",sep=""))

p <- ggplot(e,aes(x=reorder(VARIABLE, order), y=-log10(P), fill=VARIABLE)) + 
  geom_bar(stat="identity", width=.7, position = position_dodge(width = 0.6), aes(fill = ct_color)) + 
  coord_flip() + 
  geom_hline(yintercept =-log10(0.05/length(unique(e$VARIABLE))), linetype="dashed", alpha=.5) + 
  xlab("") + ylab(expression('-log'[10]*'(pvalue)')) + guides(fill="none") +
  scale_fill_identity()+
  theme(panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "darkgray", 
                                 size = .5, linetype = "solid"))

ggsave(p, 
       file=paste(mypath,"p_h2SpecificGenes.pdf", sep=""),
       width=5, height=4)
ggsave(p, 
       file=paste(mypath,"p_h2SpecificGenes.png", sep=""),
       width=5, height=4, dpi = 600)
