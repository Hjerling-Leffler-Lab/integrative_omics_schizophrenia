# 02_plot_fig17_DEG_Deciles_16CT.R

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("./7_heritability_enrichment/output/figS17/")
mydf1 <- fread("df_figS17_DEG_Deciles_16CT.tsv")

#- plot
p1a <- ggplot(mydf1 %>% 
                filter(substr(CT1,1,3)=="Exc"),
              aes(x=Decile, y=-log10(P), fill=up_down, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(CT1~up_down) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = seq(1, 10, by = 1)) +
  scale_fill_manual(values=c("#008080","darkred"),name="")
p1b <- ggplot(mydf1 %>% 
                filter(substr(CT1,1,3)=="Inh"),
              aes(x=Decile, y=-log10(P), fill=up_down, alpha=if.sig.fdr)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(CT1~up_down) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = seq(1, 10, by = 1)) +
  scale_fill_manual(values=c("#008080","darkred"))
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
p1 <- (p1a+guides(alpha=guide_legend(title="FDR sig."))) | 
  ((p1b+guides(alpha="none",fill="none")+xlab(""))/
     (p1c+guides(alpha="none",fill="none")))

#- output
ggsave(p1, 
       filename="p_figS17_DEG_Deciles_16CT.pdf",
       width = 8,height=8)

#-- end --#