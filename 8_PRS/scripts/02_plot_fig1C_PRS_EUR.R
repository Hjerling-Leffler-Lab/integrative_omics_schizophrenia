#- 02_plot_fig1C_PRS_EUR.R
#- by Shuyang Yao (shuyang.yao@ki.se)
#- 2024-10-31
#- purpose: generate figure 1C: box plot for standardized schizophrenia polygenic risk score (PRS) in samples with European ancestry.

library(ggplot2)
library(ggsignif)
library(data.table)

OUTDIR <- "../output/fig1C/"
a <- fread(paste(OUTDIR,"brn_GRS_EUR.tsv",sep="")) 
mytest <- t.test(a[a$PHENO==2,"S4_0.1"],a[a$PHENO==1,"S4_0.1"]) # t-test for standardized PRS; PRS calculated with the p-threshold method with SNP-P≤0.1

p <- ggplot(a,aes(x=Phenotype, y=S4_0.1, fill=Phenotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black",size=0.6, alpha=0.5) +
  scale_fill_manual(values=c("#f59f1c","#3d699b")) +
  ylab("Standardized schizophrenia PRS") + xlab("(European ancestry)") +
  theme_classic() +
  geom_signif(comparisons=list(c("case","control")),
              annotations = paste("p=",round(mytest$p.value,4),sep=""), size=0.3, textsize = 4) +
  guides(fill="none") +
  theme(axis.title.x = element_text(size=8))

ggsave(p, file=paste(OUTDIR,"",sep="fig1C_PRS_EUR.pdf"), height = 3.5,width = 3)

