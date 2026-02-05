# 03_plot_fig4F_corr_btw_3CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# updated: 2024-01-15 
# aim: plot figure 4F

library(data.table)
library(ggplot2)
library(patchwork)

setwd("./9_sc_eQTL/output/fig4F/")

dat3 <- fread("df_3CT0.05_vs_GTEx.tsv.gz") 
# Note on the data:
# The eQTL pairs were mapped across the three cell classes. An eQTL significant in any cell class was kept in the dataset.

#- plot aFC
dat3$aFC_Ex <- 2^(dat3$effect_size_Excitatory)
dat3$aFC_In <- 2^(dat3$effect_size_Inhibitory)
dat3$aFC_NN <- 2^(dat3$effect_size_NonNeuronal)

axislow=0
axishigh=5
p4.ei <- ggplot(dat3, 
                aes(x=aFC_Ex, y=aFC_In) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat3$aFC_Ex, dat3$aFC_In),2),sep=""))
p4.en <- ggplot(dat3, 
                aes(x=aFC_Ex, y=aFC_NN) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat3$aFC_Ex, dat3$aFC_NN),2),sep=""))
p4.in <- ggplot(dat3, 
                aes(x=aFC_In, y=aFC_NN) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat3$aFC_In, dat3$aFC_NN),2),sep=""))

p4 <- (p4.ei + guides(fill="none")) | (p4.en + guides(fill="none")) |(p4.in + guides(fill="none"))

ggsave(plot=p4,
       filename=paste(outdir,"fig4F_corr_btw_3CT.pdf",sep=""),
       width=9, height = 3)
