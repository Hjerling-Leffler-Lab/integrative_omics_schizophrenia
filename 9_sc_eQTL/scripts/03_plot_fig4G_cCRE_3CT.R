# 03_plot_fig4G_cCRE_3CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# aim: plot figure 4G

library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4G/")

dat <- fread("01_df_cCRE_around_eQTL_3CT.tsv.gz") 
# Note on the data:
# from qtltools, fdensity command, the output has three columns:
# 1    	The start position of the bin
# 2    	The end position of the bin
# 3    	The number of associations in this bin

#-  plot
p <- ggplot(dat, aes(x=(V1+V2)/2000, y=V3,color=CT_qtl)) + # get the middle point of each 1kb bin by (V1+V2)/2000.
  geom_line() +
  facet_grid(CT_qtl~CT_cCRE) +
  scale_color_manual(values=c("#00ADEE","#F05A28","black"))+
  xlab("Distance to cis-eQTL (KB)") + ylab("Nr overlaps with cCREs") +
  theme_bw()
ggsave(plot = p,filename = outplot, width = 33, height = 5)

fwrite(d,file=outdf,sep="\t", col.names = T)
system(paste("rm ",outdf,".gz",sep=""))
system(paste("gzip ",outdf,sep=""))

# cCRE cell types for plotting
mylist <- c(
  "IT_L2_3", "IT_L4_5","L5_6_NP",
  "LAMP5_1","PVALB","SST","VIP",
  "MGC","OGC","OPC")
d1 <- d %>%
  filter(CT_cCRE %in% mylist)
d1$CT_cCRE = factor(d1$CT_cCRE, levels=mylist)

p <- ggplot(d1, aes(x=(V1+V2)/2000, y=V3,color=CT_qtl)) +
  geom_line() +
  facet_grid(CT_qtl~CT_cCRE) +
  scale_color_manual(values=c("#00ADEE","#F05A28","black"))+
  xlab("Distance to cis-eQTL (per KB in ±1MB window)") + ylab("Nr overlaps with cCREs") +
  theme_bw() + 
  guides(color="none")+ theme(axis.text.x=element_blank())
ggsave(plot = p,filename = "fig4G_cCRE_3CT.pdf", width = 9, height = 3.5)


#--- end ---# 