# 02_plot_fig2G_h2TopDEG_16CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-08
# aim: plot figure 2F

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readxl)

setwd("./7_heritability_enrichment/output/fig2F/")
mydf <- fread("df_h2TopDEG_16CT.tsv")

df.color <- read_excel("./Cell_type_colors.xlsx") %>%
  distinct(`Cell type (15)`, `Cell type (15) color`) %>%
  rename(CT=`Cell type (15)`,
         ct_color=`Cell type (15) color`)
df.color$CT=gsub(" ","_",df.color$CT)

tmp <- mydf%>%filter(up_down=="down")%>%distinct(CT, if.sig.fdr) 

mydf2 <- mydf %>% 
  pivot_wider(id_cols = c(CT,CT1), names_from=up_down, values_from = P) %>%
  left_join(df.color, by="CT") %>%
  left_join(tmp, by="CT") %>%
  mutate(CT1 = gsub("(?<=\\b\\d)\\s(?=\\d\\b)", "-", CT1, perl = TRUE))

p <- ggplot(mydf2, aes(x=-log10(down), y=-log10(up), color=CT1)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, color="darkgrey") + 
  #xlim(0,3.5)+ylim(0,3.5) +
  scale_color_manual(values=mydf2$ct_color) +
  geom_text_repel(data=subset(mydf2, if.sig.fdr=="yes"), 
                  aes(label=CT1),size=4, segment.color = NA,
                  position = position_nudge_repel(x = 0.01, y = 0.005),
                  show.legend = FALSE) +
  guides(color=guide_legend(title="")) +
  xlab(expression('Top down-DEG, h'^2*'-enrichment, -log'[10]*'(P)'))+
  ylab(expression('Top up-DEG, h'^2*'-enrichment, -log'[10]*'(P)'))

#- output
ggsave(p, 
       filename="p_topDEG_16CT.pdf",
       width = 9.5,height=4)

#--- end ---#