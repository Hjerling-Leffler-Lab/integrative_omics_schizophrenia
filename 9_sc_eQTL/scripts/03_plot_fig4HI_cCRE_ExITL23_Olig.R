# 03_plot_fig4HI_cCRE_ExITL23_Olig.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# aim: plot figure 4HI

library(data.table)
library(tidyverse)
library(ggplot2)
library(readxl)

setwd("./9_sc_eQTL/output/fig4HI/")

d <- fread("01_df_cCRE_around_eQTL_15CT.tsv.gz") 
# Note on the data:
# from qtltools, fdensity command, the output has three columns:
# 1    	The start position of the bin
# 2    	The end position of the bin
# 3    	The number of associations in this bin
df.color <- read_excel("./workflow/Cell_type_colors.xlsx")

#- plot IT_L2_3
a <- d %>% 
  filter(CT_qtl2 %in% c(
    "Ex_L2-3_IT_I","Ex_L2-3_IT_II","Ex_L3-4_IT","Ex_L3-6_IT","Ex_L5-6_IT_I","Ex_L5-6_IT_II" 
  )) %>%
  filter(CT_cCRE=="IT_L2_3") %>%
  left_join(df.color%>%select(color, CT2),by=c("CT_qtl2"="CT2"))
a$CT_cCRE = factor(a$CT_cCRE, levels=unique(a$CT_cCRE))
p1 <- ggplot(a, aes(x=(V1+V2)/2000, y=V3, color=CT_qtl2)) +
  geom_line() +
  facet_grid(.~CT_qtl2) +
  scale_color_manual(values=unique(a$color),name="eQTL cell types")+
  guides(color="none") + xlab("") + ylab("") +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(plot = p1,filename = "fig4H_IT_L2_3_cCRE_around_Ex_eQTLs.pdf", width = 6, height = 1.5)

#- plot OGC
a <- d %>% 
  filter(CT_qtl2 %in% c(
    "Astrocytes","Microglial","OPC","Olig"
  )) %>%
  filter(CT_cCRE=="OGC") %>%
  left_join(df.color%>%select(color, CT2),by=c("CT_qtl2"="CT2"))
a$CT_cCRE = factor(a$CT_cCRE, levels=unique(a$CT_cCRE))
p2 <- ggplot(a, aes(x=(V1+V2)/2000, y=V3, color=CT_qtl2)) +
  geom_line() +
  facet_grid(.~CT_qtl2) +
  scale_color_manual(values=unique(a$color),name="eQTL cell types")+
  guides(color="none") + xlab("") + ylab("") +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(plot = p2,filename = "fig4I_Olig_cCRE_around_NN_eQTLs.pdf", width = 4, height = 1.5)

