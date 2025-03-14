library(data.table)
library(ggplot2)

setwd("./9_sc_eQTL/output/fig4E/")

dat <- fread("df_3CT0.05_vs_JB.tsv") 
# Note on the data: 
# The data was mapped between our eQTL to the cell type significant eQTLs from Bryois et al. 2022 (PMID=35915177) by both the SNP and the gene in an eQTL pair.
# In mapping, we required the SNPs have the same reference and alternative allele (hg19) in the two summary statisicis; 
# if they were flipped between the two summary statisitcs, we flipped the ref and alt alleles in one summary statistics, and changed the sign of the effect size accordingly.
# We mapped the corresponding cell types (Excitatory with Excitatory, Inhibitory with Inhibitory, and our NonNeuronal to Bryois Oligodendrocytes).

# plot
mycor <- cor(dat$slope, dat$effect_size_JB)
p <- ggplot(dat[sample(nrow(dat)),], # randomize the dots for better visualization
            aes(x=effect_size,y=effect_size_JB, color=CT)) +
  geom_point(alpha=.1, size=.1) + 
  xlab("Current study, effect size") + ylab("Bryois 2022, effect size") +
  theme_classic() +
  scale_color_manual(values=c("#00ADEE","#F05A28","darkgrey"))+
  guides(color = "none") +
  annotate(geom="text", x=-.5, y=1.3, 
           label=paste("corr=",round(mycor,2),sep=""),color="black",size=5)

ggsave(filename="fig4E_corr_Bryois2022_3CT.pdf",
       plot=p, width = 2.5, height = 2.5)
