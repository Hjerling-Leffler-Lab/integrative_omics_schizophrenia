# fig3B_figS10C_corr_RNAprocessing_deg.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2025-11-03
# Main plot: correlation matrix: all neuronal cell types, DEGs that overlaping any RNA-procssing terms
# supplementary plot: DEG by cell type heatmap, color by effect size

library(tidyverse)
library(data.table)
library(readxl)
library(patchwork)
library(ggplot2)
library(forcats)
library(ggtext)
library(pheatmap)

WD="/Users/shuyang.yao/Proj/SCZ_human_sc/Scripts/01_revision/99_RNAprocessing_plot/"
genemxDR <- "/Users/shuyang.yao/Proj/00GenomeMaster/"
degDR <- "/Users/shuyang.yao/Proj/SCZ_human_sc/Data/DEG_results/version202505/"

#- read in data
# read genemx
df.genemx <- fread(paste(genemxDR,"geneMatrix.tsv.gz",sep=""))

# read DEG list
df.deg <- fread(paste(degDR,"T7_DEGs_per_cell_type.csv",sep=""), quote="") %>% rename(ensgid=Gene)
colnames(df.deg)[c(1,11)] <- c("Gene_short","celltype")
df.deg$celltype <- gsub("\"","",df.deg$celltype)
df.deg$ensgid <- gsub(".*_", "",df.deg$ensgid)
df.deg$Gene_short <- gsub('^\"', '', df.deg$Gene_short)
# change celltype names to the shorter version
df.deg <- df.deg %>%
  mutate(celltype0=celltype) %>%
  mutate(celltype=case_when(
    celltype0=="Astrocytes" ~ "Astrocytes",
    celltype0=="Endothelial_and_mural_cells" ~ "Endothelial/mural",
    celltype0=="Excitatory_Layer_2_3_IT_neurons_II" ~ "Exc L2-3 IT II",
    celltype0=="Excitatory_Layer_2_3_IT_neurons_I" ~ "Exc L2-3 IT I",
    celltype0=="Excitatory_Layer_3_4_IT_neurons" ~ "Exc L3-4 IT",
    celltype0=="Excitatory_Layer_3_6_IT_neurons" ~ "Exc L3-6 IT",
    celltype0=="Excitatory_Layer_5_6_CT_and_NP_neurons" ~ "Exc L5-6 CT/NP",
    celltype0=="Excitatory_Layer_5_6_IT_neurons_II" ~ "Exc L5-6 IT II",
    celltype0=="Excitatory_Layer_5_6_IT_neurons_I" ~ "Exc L5-6 IT I",
    celltype0=="Inhibitory_LAMP5_neurons" ~ "Inh LAMP5",
    celltype0=="Inhibitory_PVALB_neurons" ~ "Inh PVALB",
    celltype0=="Inhibitory_SST_neurons" ~ "Inh SST",
    celltype0=="Inhibitory_VIP_neurons" ~ "Inh VIP",
    celltype0=="Microglial_cells" ~ "Microglial",
    celltype0=="Oligodendrocytes" ~ "Oligodendrocytes",
    celltype0=="Oligodendrocyte_progenitor_cells" ~ "OPC"
  ))
df.deg <- df.deg %>%
  select(ensgid, Gene_short, celltype, padj, negative_log10_padj, log2FoldChange)

# read RNA-processing gene list
mylist <- c("GOBP_RNA_PROCESSING", 
            "GOBP_MRNA_PROCESSING", 
            "GOBP_RNA_SPLICING"
)
pattern <- paste0("(", paste(tolower(mylist), collapse = "|"), ")$")
geneset_a <- fread(paste(genemxDR, "genesetsA.tsv.gz", sep="")) %>%
  filter(str_detect(tolower(geneset), pattern)) %>%
  filter(str_detect(tolower(geneset), str_c(tolower("gobp"), collapse = "|")))
df_mygenelist <- fread(paste(genemxDR,"genesetsB.tsv.gz",sep="")) %>%
  filter(genesetID %in% geneset_a$genesetID)

# filter for genes from the deg list 
a1 <- df_mygenelist %>% distinct(ensgid) %>%
  inner_join(df.deg, by="ensgid")

colnames(a1)
unique(a1$celltype)

# 1. get DEG list
tmp <- a1 %>% filter(padj<=0.05) %>% distinct(ensgid)

# 2. get the genes in all cell types
# Wide matrix: genes x celltypes (values = log2FC)
wide.fc <- a1 %>% filter(ensgid %in% tmp$ensgid) %>%
  select(Gene_short, celltype, log2FoldChange) %>%
  distinct() %>%                      # in case of duplicates per (gene, celltype)
  pivot_wider(names_from = celltype,
              values_from = log2FoldChange)

# 3. plot correlation of log2FC of DEGs across neuronal cell types. (main figure)
# Pairwise-comparison lets each correlation use all available genes

tmp1 <- wide.fc %>% select(starts_with("Exc"), starts_with("Inh"))
cor_mat <- cor(tmp1, tmp1, use = "pairwise.complete.obs", method = "spearman")

# Order axes: hierarchical clustering based on the correlation
hc_x <- if (nrow(cor_mat) > 1) hclust(dist(cor_mat), method = "complete") else NULL
hc_y <- if (ncol(cor_mat) > 1) hclust(dist(t(cor_mat)), method = "complete") else NULL
x_order <- rownames(cor_mat)[if (is.null(hc_x)) seq_len(nrow(cor_mat)) else hc_x$order]
y_order <- colnames(cor_mat)[if (is.null(hc_y)) seq_len(ncol(cor_mat)) else hc_y$order]

# Long format for ggplot
df_cor <- as.data.frame(as.table(cor_mat)) %>%
  mutate(
    Var1 = factor(Var1, levels = x_order),
    Var2 = factor(Var2, levels = y_order)
  )

# Plot
p1 <- ggplot(df_cor, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3) +   # << Add correlation values
  coord_fixed() +
  scale_fill_gradient2(
    limits = c(-0.02, 1),
    low = "pink", mid = "white", high = "green4",
    midpoint=0.5,
    name = "Spearman r"
  ) +
  scale_y_discrete(limits = rev(y_order)) +  # put first cluster at top
  labs(x = "", y = "",
       title = "") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  )

p1

# save data
df_cor.fc.rna.allct <- df_cor
fwrite(df_cor.fc.rna.allct, file=paste(WD,"df_02_corr_RNAprocessing_deg_allct.tsv",sep=""),
       sep = "\t", col.names = T)
# save plot
ggsave(p1, file=paste(WD,"02_corr_RNAprocessing_deg_allct.pdf",sep=""),
       width = 6, height = 5)

# 4. Heatmap of gene x cell types (all, no exclusion), based on p value (-log10padj)
# Wide matrix: genes x celltypes (values = p)
wide.fc <- as.matrix(wide.fc)
rownames(wide.fc) <- wide.fc[,1]
wide.fc <- wide.fc[,-1]
# convert to numerical
mat <- apply(wide.fc, 2, as.numeric)
rownames(mat) <- rownames(wide.fc)
colnames(mat) <- colnames(wide.fc)

# add significance label
wide.fc.long <- a1 %>% filter(ensgid %in% tmp$ensgid) %>%
  select(Gene_short, celltype, log2FoldChange, padj) %>%
  distinct() %>%
  mutate(if.sig=ifelse(padj<=0.05, "*",""))
sig_mat <- wide.fc.long %>%
  select(Gene_short, celltype, if.sig) %>%
  pivot_wider(names_from = celltype, values_from = if.sig) %>%
  as.data.frame()
# convert Gene_short to rownames
rownames(sig_mat) <- sig_mat$Gene_short
sig_mat$Gene_short <- NULL
sig_mat <- sig_mat[, colnames(mat)]
# fill NA with ""
sig_mat[is.na(sig_mat)] <- ""

# plot heatmap, color=log2FC
p2 <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  na_col = "darkgrey",
  color = colorRampPalette(c("#008080", "white", "darkred"))(100),
  main = "",
  display_numbers = sig_mat,        # add significance label
  number_color = "black",           
  fontsize_number = 10  ,
  name="log2(FC)"            
)

# save data
df_heatmap.deg.fc <- wide.fc.long
fwrite(df_heatmap.deg.fc, file=paste(WD,"df_03_headmap_rna_FC_deg_allct.tsv",sep=""),
       sep = "\t", col.names = T)
# save plot
ggsave(p2, file=paste(WD,"03_headmap_rna_FC_deg_allct.pdf",sep=""),
       width = 4.3, height = 6)
ggsave(p2, file=paste(WD,"03_headmap_rna_FC_deg_allct.png",sep=""),
       width = 4.3, height = 6, dpi=300)

#-- end --#