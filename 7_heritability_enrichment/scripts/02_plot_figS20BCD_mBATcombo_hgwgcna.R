# 02_plots_Bordeaux_mBATcombo.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-06-16
# aim: make plots of the mBATcombo results for pathways in Bordeaux modules

library(data.table)
library(tidyverse)
library(ggplot2)
library(ComplexUpset)


#--- fig 20B & D ---#
setwd("./7_heritability_enrichment/output/fig5G/")
df <- fread("df_03_pLDSC_29modules_with2moduleBlocks_shortnames.tsv")

p1 <- ggplot(df,
             aes(x=reorder(module_shortnames,-order), y= `Prop._SNPs`, fill=as.character(if.EnrichmentFDR.05)))+
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("grey","salmon"), name="FDR sig.", labels=c("no","yes")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 0.01, linetype="dashed", color="black") 

p2 <- ggplot(df,
             aes(x=reorder(module_shortnames,-order), y= -log10(Enrichment_p), fill=as.character(if.EnrichmentFDR.05)))+
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("grey","salmon"), name="FDR sig.", labels=c("no","yes")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 0.01, linetype="dashed", color="black") 

setwd("./7_heritability_enrichment/output/figS20BC/")
ggsave(p1,
       filename="FigS20B.LDSC_module_size.pdf",
       width = 4, height = 4)
ggsave(p2,
       filename="FigS20D.LDSC_enrichment_p.pdf",
       width = 4, height = 4)

#--- fig 20C ---#
setwd("./7_heritability_enrichment/output/figS20BC/")
d0 <- fread("df_figS20BC_mBATcombo.tsv")
df1 <- d0 %>% 
  filter(FDR_mBATcombo_allGenes40Modules<=0.05) %>%
  distinct(ENSGID, module_shortnames) %>%
  mutate(value=1) %>%
  pivot_wider(id_cols = ENSGID, names_from=module_shortnames, values_from = value)
df1[is.na(df1)] <- 0
module_shortnames <- colnames(df1)[-1]
p1.bmSig <- upset(df1, module_shortnames, 
                  name='Significant genes', 
                  width_ratio=0.25,
                  themes=upset_modify_themes(
                    list(
                      'overall_sizes'=theme(
                        axis.text=element_text(angle=0)
                      ))
                  )
)
ggsave(p1.bmSig,
       filename="figS20C_upset_BordeauxSigGene.pdf",
       width = 6, height = 5)

#-- end --#
