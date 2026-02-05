#- figS25_DEG_medication.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2025-11-17

library(data.table)
library(tidyverse)
library(ggplot2)
library(readxl)

# paths
path.deg <- "T7_DEGs_per_cell_type.csv" # table 5, DEG results
setwd("./Pathway_specific/figS25_medication_gene_overlap/output/")
path.med <- "01_processed_36932057_supptab3.xlsx" # supplementary table 3 from PMID 36932057.
out1p <- "figS25A_n_gene_overlap_DEG_drug.pdf"
out2p <- "figS25C_corr_per_celltype.pdf"
out3p <- "figS25B_corr_all_celltypes.pdf"

# read data
df.deg <- fread(path.deg) %>% separate(Gene, into=c("gene_name","ensgid"),sep="_",remove=F)
colnames(df.deg)[13] <- "celltype"
df.deg$celltype <- gsub("\"","",df.deg$celltype)
df.med <- read_excel(path.med)
colnames(df.med)[3] <- "ensgid"
df.med$Drug <- gsub("CLZ","Clozapine",df.med$Drug)
df.med$Drug <- gsub("HAL_hi","Haloperidol, high",df.med$Drug)
df.med$Drug <- gsub("HAL_lo","Haloperidol, low",df.med$Drug)

a <- df.deg %>% 
  filter(padj<=0.05) %>% 
  select(ensgid, gene_name, log2FoldChange, padj, celltype) %>% 
  rename(log2fc_deg=log2FoldChange, padj_deg=padj) %>%
  inner_join(df.med %>% 
               filter(is.na(ensgid)==F) %>%
               select(ensgid, logFC, P, FDR, Drug) %>%
               rename(log2fc_drug=logFC,
                      praw_drug=P,
                      padj_drug=FDR), 
             by="ensgid")
# check duplication:
tmp <- a %>% group_by(celltype, Drug, ensgid) %>% summarise(n=n()) # ENSG00000134825 mapped to 2 Macaque genes that was included in the haloperidoal-high group, the 2 genes had the same direction of logFC, keep the more significant one
a <- a %>% group_by(celltype, Drug, ensgid) %>% slice(1) %>% ungroup()
a$celltype <- gsub("_"," ",a$celltype)
# check number of overlapping genes between current DEG and the drug-induced diffrentially expressed genes.
b <- a %>% distinct(celltype, Drug, ensgid) %>% group_by(celltype, Drug) %>% summarise(n=n()) %>% ungroup()

# plot number of overlapping genes
b2 <- b %>% 
  pivot_wider(names_from=Drug, values_from=n, values_fill=0) %>%
  pivot_longer(cols = -celltype, names_to = "Drug", values_to = "n") %>%
  mutate(n_cat = ifelse(n >= 10, "≥10", "<10"))
p1 <- ggplot(b2, aes(x = Drug, y = celltype, fill = n_cat)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n), size = 3) +
  scale_fill_manual(
    values = c(
      "≥10" = "salmon",
      "<10" = "grey90"
    ),
    name = "Nr.overlapping\ngenes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Drug", y = "Cell type")
p1
ggsave(p1, file=out1p, width = 5.5,height = 3)

# keep combinations of cell-type and drug, calculate correlation of logFC and plot
b1 <- b %>% mutate(group=paste(celltype, Drug, sep="\n ")) %>% filter(n>=10)
c <- a %>% 
  mutate(group=paste(celltype, Drug, sep="\n ")) %>% 
  filter(group%in%b1$group) 

library(rstatix)
c1 <- c %>%
  group_by(group) %>%
  summarise(cor=cor(log2fc_deg, log2fc_drug, method = "spearman"),
            p   = cor.test(log2fc_deg, log2fc_drug, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(fdr=adjust_pvalue(p)) %>%
  mutate(label=paste("corr=",round(cor,2),", FDR=", signif(fdr,2),sep=""))

p2 <- ggplot(c, aes(x = log2fc_deg, y = log2fc_drug)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3, color="grey") +
  geom_vline(xintercept = 0, linewidth = 0.3, color="grey") +
  theme_bw()+
  geom_smooth(method="lm", col="black")+
  xlab("log2(Fold-change), current study DEG") + ylab("log2(Fold-change), drug") + 
  facet_wrap(~group, nrow=1) +
  geom_text(
    data = c1,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.2, vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  )
p2
ggsave(p2, file=out2p, height=2.5, width = 12)

# plot all 
d1 <- a %>%
  mutate(cor=cor(log2fc_deg, log2fc_drug, method = "spearman"),
         p   = cor.test(log2fc_deg, log2fc_drug, method = "spearman", exact = FALSE)$p.value,
         .groups = "drop"
  ) %>%
  mutate(label=paste("corr=",round(cor,2),", p=", signif(p,2),sep="")) %>%
  distinct(label)

p3 <- ggplot(a, aes(x = log2fc_deg, y = log2fc_drug)) +
  geom_point(alpha = 0.6, size = 1,
             aes(color=celltype, shape=Drug)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color="grey") +
  geom_vline(xintercept = 0, linewidth = 0.3, color="grey") +
  theme_bw()+
  geom_smooth(method="lm", col="grey40")+
  xlab("log2(Fold-change), current study DEG") + ylab("log2(Fold-change), drug") + 
  #facet_wrap(~group) +
  geom_text(
    data = d1,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  )
ggsave(p3, file=out3p, height=4, width = 7)

#--- end ---#