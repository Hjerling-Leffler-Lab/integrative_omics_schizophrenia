# 03_plot_figS14_cCRE_16CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# aim: plot figure S14

library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./9_sc_eQTL/output/figS14/")

d <- fread("01_df_cCRE_around_eQTL_15CT.tsv.gz") 
# Note on the data:
# from qtltools, fdensity command, the output has three columns:
# 1    	The start position of the bin
# 2    	The end position of the bin
# 3    	The number of associations in this bin

#- plot
p <- ggplot(d, aes(x=(V1+V2)/2000, y=V3,color=CT_qtl)) +
  geom_line() +
  facet_grid(CT_qtl2~CT_cCRE) +
  scale_color_manual(values=unique(d$color),name="eQTL cell types")+
  xlab("Distance to cis-eQTL (KB)") + ylab("Nr overlaps with cCREs") +
  theme_bw() + theme(axis.text.x=element_blank())
ggsave(plot = p,filename = "figS14_cCRE_16CT.pdf", width = 30, height = 16)

#--- end ---# 