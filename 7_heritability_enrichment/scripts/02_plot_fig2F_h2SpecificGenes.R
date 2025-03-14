# 02_plot_fig2f.R

library(data.table)
library(tidyverse)
library(ggplot2)

# path
mypath <- "./7_heritability_enrichment/output/fig2F/"

e <- fread(paste(mypath,"df_fig2F_h2SpecificGenes.tsv",sep=""))

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
       file=paste(mypath,"p_fig2F_h2SpecificGenes.pdf", sep=""),
       width=5, height=4)
ggsave(p, 
       file=paste(mypath,"p_fig2F_h2SpecificGenes.png", sep=""),
       width=5, height=4, dpi = 600)
