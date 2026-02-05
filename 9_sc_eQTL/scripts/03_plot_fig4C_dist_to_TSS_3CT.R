# 03_plot_fig4C_dist_to_TSS_3CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-04-30 
# aim: plot figure 4C

library(data.table)
library(tidyverse)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4C/")
datapath <- "./9_sc_eQTL/output/3CT/"

ex <- fread(paste(datapath,"Excitatory_all.chr.eQTL_combined_qvalue0.05.tsv.gz",sep="")) %>%
  select(phe_id, var_id, dist_phe_var) %>%
  mutate(cell_type="Excitatory")
inh <- fread(paste(datapath,"Inhibitory_all.chr.eQTL_combined_qvalue0.05.tsv.gz",sep=""))  %>%
  select(phe_id, var_id, dist_phe_var) %>%
  mutate(cell_type="Inhibitory")
nn <- fread(paste(datapath,"NonNeuronal_all.chr.eQTL_combined_qvalue0.05.tsv.gz",sep=""))  %>%
  select(phe_id, var_id, dist_phe_var) %>%
  mutate(cell_type="NonNeuronal")

dat <- rbind(ex, inh, nn) %>%
  mutate(dist_phe_var2=dist_phe_var/1000) # define distance in KB.

p <- ggplot(dat, aes(x=dist_phe_var2,color=cell_type, fill=cell_type)) +
  geom_density(alpha=.05) +
  theme(panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = "darkgray", 
                                 size = .5, linetype = "solid"),
        axis.ticks.y = element_blank()) +
  xlab("Distance to TSS (Kb)") + ylab("Density") +
  scale_color_manual(values=c("#00ADEE","#F05A28","#808080")) +
  scale_fill_manual(values=c("#00ADEE","#F05A28","#808080")) +
  guides(fill="none")
p

ggsave(p, file="fig4C_dist_to_TSS_3CT.pdf",
       width = 5, height=2)
