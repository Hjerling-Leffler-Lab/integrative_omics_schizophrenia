# 03_plot_fig4E_corr_Bryois.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# updated: 2024-01-15 
# aim: plot figure 4E

library(data.table)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4E/")

dat <- fread("df_3CT0.05_vs_JB.tsv") 
# Note on the data: 
# The data was mapped between our eQTL to the cell type significant eQTLs from Bryois et al. 2022 (PMID=35915177) by both the SNP and the gene in an eQTL pair.
# In mapping, we required the SNPs have the same reference and alternative allele (hg19) in the two summary statisicis; 
# if they were flipped between the two summary statistics, we flipped the ref and alt alleles in one summary statistics, and changed the sign of the effect size accordingly.
# We mapped the corresponding cell types (Excitatory with Excitatory, Inhibitory with Inhibitory, and our NonNeuronal to Bryois Oligodendrocytes).

# derive allele fold change
dat$aFC <- 2^(dat$slope)
dat$aFC_JB <- 2^(dat$effect_size_JB)

# plot
mycor <- cor(dat$aFC, dat$aFC_JB) #0.8465326
p <- ggplot(dat[sample(nrow(dat)),], 
             aes(x=aFC,y=aFC_JB, color=CT)) +
  geom_point(alpha=.1, size=.1) + 
  xlab("Current study, aFC") + ylab("Bryois 2022, aFC") +
  theme_classic() +
  scale_color_manual(values=c("#00ADEE","#F05A28","darkgrey"))+
  guides(color = "none") +
  annotate(geom="text", x=1.2, y=3.2, 
           label=paste("corr=",round(mycor,2),sep=""),color="black",size=5)

ggsave(filename="fig4E_corr_Bryois2022_3CT.pdf",
       plot=p, width = 2.5, height = 2.5)
