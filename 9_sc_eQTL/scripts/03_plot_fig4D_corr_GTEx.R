# 03_plot_fig4D_corr_GTEx.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2024-10-24 
# updated: 2024-01-15 
# aim: plot figure 4D

library(data.table)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4D/")

dat <- fread("df_3CT0.05_vs_GTEx.tsv") 
# Note on the data: 
# The data was mapped between our eQTL to GTEx BrainPFC (BA9) significant eQTLs by both the SNP and the gene in an eQTL pair.
# In mapping, we required the SNPs have the same reference and alternative allele (hg19) in the two summary statistics; 
# if they were flipped between the two summary statistics, we flipped the ref and alt alleles in one summary statistics, and changed the sign of the effect size accordingly.

# convert beta to allele frequency
dat$aFC <- 2^(dat$effect_size)
dat$aFC_GTEx <- 2^(dat$effect_size_GTEx)

# plot
mycor <- cor(dat$aFC, dat$aFC_GTEx) #0.8379113
p <- ggplot(dat[sample(nrow(dat)),], 
             aes(x=aFC,y=aFC_GTEx, color=CT)) +
  geom_point(alpha=.1, size=.1) + 
  xlab("Current study, aFC") + ylab("GTEx-v8 BA9, aFC") +
  theme_classic() +
  scale_color_manual(values=c("#00ADEE","#F05A28","darkgrey"))+
  guides(color = "none") +
  annotate(geom="text", x=1, y=3, 
           label=paste("corr=",round(mycor,2),sep=""),color="black",size=5)

ggsave(filename="fig4D_corr_GTEx_3CT.pdf",
       plot=p, width = 2.5, height = 2.5)
