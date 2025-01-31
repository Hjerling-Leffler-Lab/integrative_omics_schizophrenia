#- 02_plot_fig1B_genotypePCs.R
#- Shuyang Yao, 20210525
#- goal: plot genotyping PCs (fig1B)
#- last updated: 20220613

library(tidyverse)
library(data.table)
library(ggplot2)
library(readxl)

#- path
OUTDIR <- "./8_PRS/output/fig1B/"

#- sample info (supplementary table 1)
tab1 <- read_excel("./workflow/Supplementary_Table1.xlsx") 
tab1rm <- tab1 %>% filter(snRNAseq_excluded_after_QC=="yes")

#- read in PC data
dmerge <- fread(paste(OUTDIR,"brn_pca_noref.menv.mds",sep="")) 

#- merge for sample information
dmerge <- dmerge %>%
  mutate(ancestry=ifelse((C1>-0.2 & C2>-0.1),"EUR","non-EUR")) %>% # define EUR ancestry
  filter(!IID%in%tab1rm$donor_ID_genotype) %>% # remove the two samples that did not path snRNA-seq QC
  left_join(tab1%>%select(donor_ID_genotype, scz_status),by=c("IID"="donor_ID_genotype"))

# gerenate dataset with sample ancestry info
fwrite(dmerge %>% select(IID, ancestry), file="./workflow/sample_brn_ancestry.tsv", sep="\t", col.names = T)

#- plot
p <- ggplot(dmerge, aes(C1, C2, colour=scz_status, shape=ancestry)) + 
  geom_point(alpha=.9, size=2) + 
  scale_color_manual(values=c("#f59f1c","#3d699b")) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent"), 
        panel.grid = element_blank(),
        axis.text = element_blank()
  ) +
  labs(x="PC1", y="PC2") 

ggsave(p, file=paste(OUTDIR,"fig1B_genotypePCs.pdf",sep=""), height = 3,width = 4)

#--- ends ---#