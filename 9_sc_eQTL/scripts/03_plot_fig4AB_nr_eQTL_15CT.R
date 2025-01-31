library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4AB/")
mydf <- fread("01_df_eQTL_numbers_16CT.tsv")

#- fig 4A: number of eGenes
p1 <- ggplot(mydf, aes(x=CT, y=n_eQTL, fill=CT)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=rev(mydf$color)) +
  coord_flip() +
  theme_classic() +
  ylab("Nr eGene")

#- fig 4B: number of cells and number of eGenes
p2 <- ggplot(mydf, aes(x=n_cell, y=n_eGene, color=CT)) +
  geom_point() +
  scale_color_manual(values=rev(b2$color)) +
  coord_flip() +
  theme_classic()

#- output
ggsave(p1+guides(fill="none"),filename="fig4A_bar_NeGene.pdf",width=7,height = 5)
ggsave(p2+guides(color="none"),filename="fig4B_point_NeGene_Ncell.pdf",width=3,height = 3)

#-- end --#