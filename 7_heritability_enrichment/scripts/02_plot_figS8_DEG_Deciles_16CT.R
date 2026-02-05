# 02_plot_figS8_DEG_Deciles_16CT.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-08
# aim: plot figure S8

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("./7_heritability_enrichment/output/figS8/")
mydf1 <- fread("df_DEG_Deciles_16CT.tsv")

mydf1 <- mydf1 %>%
  mutate(alpha_val = ifelse(if.sig.fdr=="yes", 1, 0.2))  # assuming if.sig.fdr is logical

#- Excitatory cell types
p1a <- ggplot(mydf1 %>% 
                filter(substr(CT1,1,3)=="Exc"),
              aes(x=Decile, y=-log10(P), fill=up_down, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(CT1~up_down) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = seq(1, 10, by = 1)) +
  scale_fill_manual(values=c("#008080","darkred"),name="")
#- Inhibitory cell types
p1b <- ggplot(mydf1 %>% 
                filter(substr(CT1,1,3)=="Inh"),
              aes(x=Decile, y=-log10(P), fill=up_down, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(CT1~up_down) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = seq(1, 10, by = 1)) +
  scale_fill_manual(values=c("#008080","darkred"))
#- Nonneuronal cell types
p1c <- ggplot(mydf1 %>% 
                filter(CT1 %in%
                         c("Astrocytes","Endothelial/mural",
                           "Microglial cells","OPC","Oligodendrocytes")) %>%
                mutate(CT2=case_when(
                  CT1=="Astrocytes" ~ "Astr",
                  CT1=="Endothelial/mural" ~ "Endo/m",
                  CT1=="Microglial cells" ~ "Micr",
                  CT1=="Oligodendrocytes" ~ "Olig",
                  CT1=="OPC" ~ "OPC"
                )),
              aes(x=Decile, y=-log10(P), fill=up_down, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(CT2~up_down) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = seq(1, 10, by = 1)) +
  scale_fill_manual(values=c("#008080","darkred"))

# Combine plots
# p1 <- p1a | (p1b / p1c)
right_panel <- plot_grid(p1b + theme(legend.position = "none"), 
                         p1c + theme(legend.position = "none"), 
                         ncol = 1, align = "v")

# Extract legend from p1a
legend <- cowplot::get_legend(p1a + guides(alpha=guide_legend(title="FDR sig.")))

# Remove legend from p1a now (to avoid duplication)
p1a_nolegend <- p1a + theme(legend.position = "none")

# Combine left, right, and legend horizontally
final_plot <- plot_grid(p1a_nolegend, right_panel, legend, 
                        ncol = 3, 
                        rel_widths = c(2, 2, 0.4),
                        align = "h")
                        
#- output
ggsave(final_plot, 
       filename="p_figS8_DEG_Deciles_16CT.pdf",
       width = 8,height=8)

#-- end --#