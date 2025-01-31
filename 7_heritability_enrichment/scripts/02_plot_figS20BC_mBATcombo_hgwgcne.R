# 02_plots_Bordeaux_mBATcombo.R
# aim: make plots of the mBATcombo results for pathways in Bordeaux modules

library(data.table)
library(tidyverse)
library(ggplot2)
library(ComplexUpset)

setwd("./7_heritability_enrichment/output/figS20BC/")
d0 <- fread("df_figS20BC_mBATcombo.tsv")

# fig 19b, upset plot
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
       filename="figS20B_upset_BordeauxSigGene.pdf",
       width = 6, height = 5)

# fig 5g, dot plot of pathway proportion of mBATcombo sig genes
library(scatterpie)
df2 <- d0 %>%
  distinct(go_category, module_shortnames, ENSGID, if.fdr.sig_allGenes40Modules) %>%
  group_by(go_category, module_shortnames) %>%
  summarise(n=n(), sig=sum(if.fdr.sig_allGenes40Modules)) %>%
  ungroup() %>%
  mutate(not_sig=n-sig)
tmp1 <- df2 %>% distinct(go_category) %>% 
  arrange(go_category) %>% mutate(go_category2=row_number()*2)
tmp2 <- df2 %>% distinct(module_shortnames) %>%
  arrange(module_shortnames) %>% mutate(module_shortnames2=row_number()*2)
df2a <- df2 %>% 
  left_join(tmp1, by="go_category") %>% 
  left_join(tmp2, by="module_shortnames") %>%
  mutate(pct=sig/n)

p1 <- ggplot()+ 
  geom_scatterpie(data = df2a, 
                  aes(y=go_category2, x=module_shortnames2, r=log10(n)*0.32), 
                  alpha=0.9, 
                  cols = c("sig","not_sig"), color=NA) + 
  theme_classic() +
  scale_y_continuous(labels=tmp1$go_category, 
                     breaks=seq(2,2*nrow(tmp1),2), limits=c(0,nrow(tmp1)*2+0.8)) +
  scale_x_continuous(labels=tmp2$module_shortnames, 
                     breaks=seq(2,2*nrow(tmp2),2), limits=c(0,nrow(tmp2)*2+2)) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("") +
  guides(fill=guide_legend(title="mBATcombo")) +
  geom_scatterpie_legend(breaks = seq(log10(10)*0.32,
                                      log10(1000)*0.32,
                                      length=3) , 
                         x = 2, y = 2, #n=3, 
                         labeller = function(x) ceiling(10^(x/0.32)),
                         size=3) +
  coord_equal() 

ggsave(p1, filename="figS20C.dot_pie_mbatcombo_go_bordeaux.pdf",
       width=7,height = 12)

#-- end --#
