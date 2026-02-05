# fig3C_mitochondria_plot.R
# author: Shuyang Yao (shuyang.yao@ki.se)
# date: 2025-11-03
# To plot DEG results in nuclear genes related to mitonchondria subunites

library(tidyverse)
library(data.table)
library(readxl)
library(patchwork)
library(ggplot2)
library(forcats)
library(ggtext)

#--- 1. read in data ---# 
# genemx
df.genemx <- fread("geneMatrix.tsv.gz")

# DEG results
df.deg <- fread("T7_DEGs_per_cell_type.csv", quote="") %>% rename(ensgid=Gene)
# fix format
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

# proteomics results
proteinpath <- "diagnosis_vs_layers.xlsx"
df.pro <- read_excel(proteinpath) %>%
  left_join((df.genemx %>% select(ensgid, gene_name)),by=c("hgnc_symbol"="gene_name"))

# mitochondria nuclear gene list
# downloaded at https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls
df.mito0 <- read_excel("Human.MitoCarta3.0.xls", sheet="A Human MitoCarta3.0")

# filter out genes from the mito gene list that were not included in either DEG or protein list
tmp <- rbind(df.deg%>%distinct(ensgid), df.pro%>%distinct(ensgid))
df.mito <- df.mito0 %>%
  separate_rows(EnsemblGeneID_mapping_version_20200130, sep = "\\|") %>%
  filter(EnsemblGeneID_mapping_version_20200130%in%tmp$ensgid) # filter away those not tested in DEG or proteiomics
rm(tmp)

#--- 2. get pathway categoreis from mitochondria list ---# 

# 7 major functional pathway categories:
# (i) mitochondrial central dogma, 
# (ii) protein import, sorting and homeostasis, 
# (iii) oxidative phosphorylation (OXPHOS), 
# (iv) metabolism, 
# (v) small molecule transport, 
# (vi) mitochondrial dynamics and surveillance, 
# (vii) signaling

a0 <- df.mito %>% 
  mutate(pathway_cat=str_extract(MitoCarta3.0_MitoPathways, "^[^>|]+")) %>%
  mutate(pathway_cat=str_trim(pathway_cat))

# get gene list for oxidative phosphorylation components
a3 <- a0 %>% 
  filter(pathway_cat=="OXPHOS") %>%
  mutate(OXPHOS_complex=str_extract(MitoCarta3.0_MitoPathways, "(?<=>)[^>|]+(?=[>|])"))

#--- 3. Merge DEG and proteomics results to the mitochondria gene lists ---#
b3.1 <- a3 %>%
  select(Symbol, EnsemblGeneID_mapping_version_20200130, OXPHOS_complex) %>%
  left_join(
    df.deg %>% 
      select(ensgid, celltype, padj, negative_log10_padj, log2FoldChange),
    by=c("EnsemblGeneID_mapping_version_20200130"="ensgid")
  ) 
b3.2 <- a3 %>%
  select(Symbol, EnsemblGeneID_mapping_version_20200130, OXPHOS_complex) %>%
  left_join(
    df.pro %>% 
      select(ensgid, qvalue_diagnosis, log2fc_diagnosis) %>%
      mutate(negative_log10_padj=-log10(qvalue_diagnosis),
             celltype="Proteomics") %>%
      rename(padj=qvalue_diagnosis, log2FoldChange=log2fc_diagnosis) %>%
      select(ensgid, celltype, padj, negative_log10_padj, log2FoldChange),
    by=c("EnsemblGeneID_mapping_version_20200130"="ensgid")
  ) 

b3 <- rbind(b3.1,b3.2) %>% 
  rename(Gene=EnsemblGeneID_mapping_version_20200130,
         Gene_short=Symbol) %>%
  filter(is.na(celltype)==F,
         OXPHOS_complex!=" OXPHOS assembly factors ") # only the five complex

#--- 4. Plot heatmap per complex. ---#
df <- b3 

# 1. Cluster genes (Gene_short) within each OXPHOS_complex using -log10(padj)
clustered_gene_orders <- df %>%
  group_by(OXPHOS_complex) %>%
  summarise(gene_order = list({
    wide <- pivot_wider(cur_data(),
                        names_from = celltype,
                        values_from = negative_log10_padj,
                        values_fill = 0)
    mat <- as.matrix(wide[, -c(1:2)])  # drop Gene and Gene_short
    rownames(mat) <- wide$Gene_short
    hc <- hclust(dist(mat))
    rownames(mat)[hc$order]
  }), .groups = "drop")

# 2. Join clustered gene order into df
df_ordered <- df %>%
  left_join(clustered_gene_orders, by = "OXPHOS_complex") %>%
  group_by(OXPHOS_complex) %>%
  mutate(Gene_short = factor(Gene_short, levels = unique(gene_order[[1]]))) %>%
  ungroup()

# 3. Set aesthetics: color, alpha, significance
df_ordered <- df_ordered %>%
  mutate(
    color = ifelse(log2FoldChange < 0, "#008080", "darkred"),
    alpha_val = negative_log10_padj,
    significant = padj <= 0.05,
    sig_status = ifelse(significant, "Significant", "Not significant")
  )

# 4. Reorder cell types: Exc → Inh → Other
celltype_order <- df_ordered %>%
  distinct(celltype) %>%
  mutate(group = case_when(
    grepl("^Exc", celltype) ~ "1_Excitatory",
    grepl("^Inh", celltype) ~ "2_Inhibitory",
    grepl("^Prot", celltype) ~ "0_Proteomics",
    TRUE ~ "3_NonNeuronal"
  )) %>%
  arrange(group, celltype) %>%
  pull(celltype)

df_ordered$celltype <- factor(df_ordered$celltype, levels = rev(celltype_order))

# 5. color label gene names of those significant in both proteomics and DEG
color_gene_name <- df %>%
  mutate(test_type=ifelse(celltype=="Proteomics","pro","deg"),
         sig=ifelse(padj<=0.05,1,0)) %>%
  group_by(Gene_short,test_type) %>%
  summarise(sigmax=max(sig)) %>%
  ungroup() %>%
  group_by(Gene_short) %>%
  summarise(sum_sigmax=sum(sigmax)) %>%
  filter(sum_sigmax==2)

# Mark special genes for coloring
df_ordered <- df_ordered %>%
  mutate(Gene_short_label = ifelse(
    Gene_short %in% color_gene_name$Gene_short,
    paste0("<span style='color:#008080;'>", Gene_short, "</span>"),
    as.character(Gene_short)
  ))
# Make it a factor with the same order as Gene_short
df_ordered$Gene_short_label <- factor(df_ordered$Gene_short_label,
                                      levels = unique(df_ordered$Gene_short_label))

# 6. Dummy row for linetype legend (after Gene_short_label is defined)
example_gene_label <- df_ordered$Gene_short_label[which(!is.na(df_ordered$Gene_short_label))][1]
example_cell <- df_ordered$celltype[which(!is.na(df_ordered$celltype))][1]
example_complex <- df_ordered$OXPHOS_complex[1]

df_legend <- data.frame(
  sig_status = c("Significant", "Not significant"),
  Gene_short_label = factor(rep(example_gene_label, 2), levels = levels(df_ordered$Gene_short_label)),
  celltype = factor(rep(example_cell, 2), levels = levels(df_ordered$celltype)),
  OXPHOS_complex = example_complex
)

# 7. Plot rotated heatmap: x = Gene, y = Cell Type
p <- ggplot(df_ordered, aes(x = Gene_short_label, y = celltype)) +
  # Main heatmap tiles
  geom_tile(aes(fill = color, alpha = alpha_val), color = "white") +
  
  # Outline for significant DEGs (black frame)
  geom_tile(aes(linetype = sig_status), fill = NA, color = "black", linewidth = 0.5, show.legend = FALSE) +
  
  # Dummy tile for legend key
  geom_tile(data = df_legend,
            aes(x = Gene_short_label, y = celltype, linetype = sig_status),
            fill = NA, color = "black", linewidth = 0.6,
            inherit.aes = FALSE, show.legend = TRUE) +
  
  # Scales
  scale_fill_identity() +
  scale_alpha(range = c(0.2, 1), name = "-log10 adjusted p-value") +
  
  scale_linetype_manual(
    name = "DEG, proteomics significant",
    values = c("Significant" = "solid", "Not significant" = "blank")
  ) +
  scale_y_discrete(limits = levels(df_ordered$celltype)) +
  
  # Facet by mitochondrial complex
  facet_grid(cols = vars(OXPHOS_complex), scales = "free_x", space = "free_x") +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),  # colored gene labels
    axis.text.y = element_text(size = 8),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.4, "lines"),
    strip.text.x = element_text(angle = 0, face = "bold"),
    strip.placement = "outside",
    legend.position = "right"
  ) +
  
  # Titles
  labs(
    title = "DEG, proteomics of genes in the mitochondria complexes",
    x = "Gene", y = "Cell Type"
  )

# output
fwrite(df, file="01_df_OXPHOS_complex.tsv", sep="\t", col.names = T)
ggsave(p, file="fig3C_mitochondria_plot.pdf", height = 3.2, width = 19)

#--- end ---#