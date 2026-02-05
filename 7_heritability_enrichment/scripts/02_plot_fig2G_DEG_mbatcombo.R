# 02_plot_fig2G_DEG_mbatcombo.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2025-06-25
# aim: plot figure 2G

library(data.table)
library(readxl)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(pheatmap)

setwd("./7_heritability_enrichment/output/fig2G/")

#- read GWAS gene lists (mbatcombo)
df1.mbatcombo <- fread("full_mbatcombo_scz2022_eur.tsv") %>%
  adjust_pvalue(p.col="P_mBATcombo", output.col = "FDR", method = "fdr") %>%
  adjust_pvalue(p.col="P_mBATcombo", output.col = "P_mBATcombo_bonf", method = "bonferroni")
df1a <- df1.mbatcombo %>% filter(P_mBATcombo_bonf<=0.05)

df.deg <- fread("T7_DEGs_per_cell_type.csv", quote="")
colnames(df.deg)[c(1,11)] <- c("Gene_short","celltype")
df.deg$celltype <- gsub("\"","",df.deg$celltype)
df.deg$Gene <- gsub(".*_", "",df.deg$Gene)
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
  select(Gene, Gene_short, celltype, padj, negative_log10_padj, log2FoldChange)

#- merge
df <- df1a %>% select(Gene) %>%
  inner_join(df.deg, by="Gene")

# 1. keep genes with any significant DEG results to plot heatmap  
signif_genes <- df %>%
  filter(padj <= 0.05) %>%
  pull(Gene) %>%
  unique()

df_filt <- df %>%
  filter(Gene %in% signif_genes)

# 2. Pivot to wide format for clustering
mat_wide <- df_filt %>%
  select(Gene_short, celltype, log2FoldChange) %>%
  pivot_wider(names_from = celltype, values_from = log2FoldChange, values_fill = 0)

log2fc_matrix <- as.matrix(mat_wide[,-1])
rownames(log2fc_matrix) <- mat_wide$Gene_short

# 3. Hierarchical clustering
gene_clust <- hclust(dist(log2fc_matrix), method = "complete")
gene_order <- rownames(log2fc_matrix)[gene_clust$order]

# 4. Reorder based on clustering
df_filt$Gene_short <- factor(df_filt$Gene_short, levels = gene_order)

# 5. specify color
df_filt <- df_filt %>%
  mutate(
    color = ifelse(log2FoldChange < 0, "#008080", "darkred"),
    alpha_val = negative_log10_padj,
    significant = padj <= 0.05
  )

# 6. Extract significant data for outline tiles
df_signif <- df_filt %>% filter(significant)

# 7. # Add linetype category
df_filt <- df_filt %>%
  mutate(
    sig_status = ifelse(significant, "Significant", "Not significant")
  )

# 8. reorder cell types.
df_filt <- df_filt %>%
  mutate(
    cell_group = case_when(
      grepl("^Exc", celltype) ~ "Excitatory",
      grepl("^Inh", celltype) ~ "Inhibitory",
      TRUE ~ "Other"
    )
  )
celltype_order <- df_filt %>%
  distinct(celltype, cell_group) %>%
  arrange(factor(cell_group, levels = c("Excitatory", "Inhibitory", "Other")), celltype) %>%
  pull(celltype)

# Apply ordered factor to both df_filt and df_signif
df_filt$celltype <- factor(df_filt$celltype, levels = celltype_order)
df_signif <- df_filt %>% filter(significant)

# 9. Plot
p <- ggplot(df_filt, aes(y = Gene_short, x = celltype)) +
  geom_tile(aes(fill = color, alpha = alpha_val), color = "white") +
  scale_fill_identity() +
  scale_alpha(range = c(0.2, 1), name = "-log10(DEG adjusted P)") +
  geom_tile(aes(linetype = sig_status), fill = NA, color = "black", linewidth = 0.5) +
  scale_linetype_manual(
    name = "DEG significant",
    values = c("Significant" = "solid", "Not significant" = "blank")
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title = "Expression of GWAS genes with significant DEG",
    x = NULL,
    y = NULL
  ) 

#- output
fwrite(df_filt, "df_DEG_GWAS.tsv", sep="\t", col.names = T)
ggsave(p, 
       filename="fig_DEG_GWAS.pdf",
       height=7, width = 5)

#-- end --#