## about: plotting functions for running hdWGCNA
## author: Lisa Bast
## date: 29.06.2023
## version: 0.0.1

### Functions for plotting 
visualize_meta_cells <- function(seurat_obj, settings, opt_data){
  seurat_obj <- ScaleMetacells(seurat_obj, 
                               features=VariableFeatures(seurat_obj))
  if (settings$dim_reduction == 'pca'){
    seurat_obj <- RunPCAMetacells(seurat_obj, 
                                  features=VariableFeatures(seurat_obj))
  } else if (settings$dim_reduction == 'harmony'){
    seurat_obj <- RunHarmonyMetacells(seurat_obj, 
                                      group.by.vars=opt_data$sample_variable)
  }
  seurat_obj <- RunUMAPMetacells(seurat_obj, 
                                 reduction=settings$dim_reduction, 
                                 dims=1:15)
  
  p1 <- DimPlotMetacells(seurat_obj, group.by=settings$level,shuffle=TRUE) + 
    umap_theme() + 
    ggtitle("Cell Type") #+
  #theme(legend.position="bottom")#, reduction=settings$dim_reduction
  
  p2 <- DimPlotMetacells(seurat_obj, group.by=opt_data$sample_variable, shuffle=TRUE) + 
    umap_theme() + 
    theme(legend.text = element_text(size=4),
          legend.key.size = unit(4, 'mm')) +
    guides(color = guide_legend(override.aes = list(size=4), ncol=2)) +
    ggtitle("Sample") #+ 
  #theme(legend.position="bottom")
  
  p3_title <- "Disease"
  
  p3 <- DimPlotMetacells(seurat_obj, group.by=opt_data$group, shuffle=TRUE) + 
    umap_theme() + 
    ggtitle(p3_title) #+ 
  #theme(legend.position="bottom")
  
  #p = plot_grid(p1,p2,p3,labels=c("A","B","C"), ncol=1, nrow=3)
  
  if (settings$wgcna_name == 'scz_metacell_test'){
    results_path <- opt_data$results_path_long
  } else{
    results_path <- opt_data$results_path
  }

  add_str <- paste0("_fraction_genes",str_replace(as.character(settings$fraction_of_cells_expr_gene),".",""))
  # all three together in one pdf
  graphic_filename <- paste0(results_path,"metacells.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 15) # The height of the plot in inches
  print(p1 / p2 / p3)
  dev.off()
}

visualize_modules <- function(seurat_obj,results_path){
  p <- PlotKMEs(seurat_obj, ncol=5, text_size =1)
  
  graphic_filename <- paste0(results_path,"modules_KMEs.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p)
  dev.off()
}

plot_module_network_with_genes_highlighted <- function(TOM,modules,gene_lists,gene_list_name,ct,module_name,results_path,highlight_colors,other_color){
  
  # get module colors for plotting 
  mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
  mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)
  
  # set up the graph object with igraph & tidygraph
  graph <- TOM %>% 
    igraph::graph_from_adjacency_matrix(mode='undirected', weighted=TRUE) %>% 
    tidygraph::as_tbl_graph(directed=FALSE) %>% 
    tidygraph::activate(nodes) 
  
  cols_f <- function(genes){
    color <- replicate(length(genes),other_color)
    i=0
    for (g in genes){
      i<-i+1
      if (gene_list_name=="all"){
        counter=0
        if (g %in% gene_lists$common_variants){
          color[i]=highlight_colors[1]
          counter = counter+1
        } 
        if (g %in% gene_lists$WES_genes){
          color[i]=highlight_colors[2]
          counter = counter+1
        }
        if (g %in% gene_lists$mBAT_genes){
          color[i]=highlight_colors[3]
          counter = counter+1
        }
        if (g %in% gene_lists$DEP_genes){
          color[i]=highlight_colors[4]
          counter = counter+1
        }
        if (counter>1){
          color[i]=highlight_colors[5]
          counter = counter+1
        }
      } else if(gene_list_name=="common_variants"){
        if (g %in% gene_lists$common_variants){
          color[i]=highlight_colors[1]
        } 
      } else if(gene_list_name=="WES_genes"){
        if (g %in% gene_lists$WES_genes){
          color[i]=highlight_colors[2]
        }
      } else if(gene_list_name=="mBAT_genes"){
        if (g %in% gene_lists$mBAT_genes){
          color[i]=highlight_colors[3]
        }
      } else if(gene_list_name=="DEP_genes"){
        if (g %in% gene_lists$DEP_genes){
          color[i]=highlight_colors[4]
        }
      }
    }
    return(color)
  }
  
  # make the plot with ggraph
  p <- ggraph(graph) + 
    geom_edge_link(color='grey', alpha=0.2) + 
    geom_node_point(color = cols_f(names(graph[1]))) +
    #geom_node_point(color='black') +#  geom_node_point(aes(colour = depth)) +
    geom_node_label(aes(label=name), repel=TRUE, max.overlaps=Inf, fontface='italic',color = cols_f(names(graph[1]))) + 
    ggtitle(paste0(ct," ",module_name," module"))
  
  graphic_filename <- paste0(results_path,"Module_Network_Plot/",ct,"_",module_name,"_",gene_list_name,".pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 5, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p)
  dev.off()
}


visualize_dendrogram <- function(seurat_obj,results_path,ct){
  #browser()
  #why is pdf corrupted? --Y cannot be stored in a variable like this
  #p <- PlotDendrogram(seurat_obj, main=paste0(ct,' hdWGCNA Dendrogram'))
  graphic_filename <- paste0(results_path,"dendrogram.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  PlotDendrogram(seurat_obj, main=paste0(ct,' hdWGCNA Dendrogram'))
  dev.off()
}

visualize_soft_powers_results <- function(seurat_obj,results_path){
  # plot the soft powers results:
  plot_list <- PlotSoftPowers(seurat_obj)
  
  # assemble with patchwork
  p <- wrap_plots(plot_list, ncol=2)
  
  graphic_filename <- paste0(results_path,"soft_powers_results.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p)
  dev.off()
}

visualize_module_features <- function(seurat_obj,results_path, features_to_visualize){
  if (features_to_visualize=='hMEs'){
    ord <- TRUE
  } else{
    ord <- "shuffle"
  }
  # make a featureplot of hMEs for each module
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features=features_to_visualize, # plot the hMEs
    order=ord, # order so the points with highest hMEs are on top
    point_size = 0.1,
  )
  # stitch together with patchwork
  p <- wrap_plots(plot_list, ncol=6)
  
  graphic_filename <- paste0(results_path,"module_features_",features_to_visualize,".pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p)
  dev.off()
}

visualize_module_correlations <- function(seurat_obj,results_path,ct){
  #browser()
  par <- c("hMEs","MEs","scores") #"average"
  for (i in seq(1,length(par))){
    p <- par[i]
    graphic_filename <- paste0(results_path,"module_correlations_",p,".pdf")
    pdf(file = graphic_filename,   # The directory you want to save the file in
        width = 12, # The width of the plot in inches
        height = 6) # The height of the plot in inches
    print(ModuleCorrelogram_v2(seurat_obj,features=p,wgcna_name = ct)) #"MEs","scores","average")
    dev.off()
  }
  
  #plot MEs vs hMEs:
  hMEs <- GetMEs(seurat_obj, harmonized=TRUE, wgcna_name=ct)
  graphic_filename <- paste0(results_path,"module_correlations_MEs_vs_hMEs.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(ModuleCorrelogram_v2(seurat_obj,MEs2=hMEs,features="MEs",wgcna_name = ct)) #"MEs","scores","average")
  dev.off()
}

visualize_hMEs_as_dotplot <- function(seurat_obj,mods,opt_data){
  #browser()
  #does it work with the disease status?
  # plot with Seurat's DotPlot function
  p1 <- DotPlot(seurat_obj, features=mods, group.by = opt_data$group)
  p2 <- DotPlot(seurat_obj, features=mods, group.by = opt_data$level)
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p1 <- p1 +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  p2 <- p2 +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  
  # save figure
  graphic_filename <- paste0(opt_data$results_path_ct_norm,"dot_plot_hMEs.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 20, # The width of the plot in inches
      height = 7) # The height of the plot in inches
  grid.arrange(p1,p2,nrow=1, widths = c(1,2))
  #print(p1 | p2)
  dev.off()
  
  return(seurat_obj)
}

visualize_hMEs_as_violins <- function(seurat_obj,level, results_path_ct, module_name){
  
  # Plot INH-M4 hME using Seurat VlnPlot function
  p <- VlnPlot(
    seurat_obj,
    features = module_name,
    group.by = level,
    pt.size = 0 # don't show actual data points
  )
  
  # add box-and-whisker plots on top:
  p <- p + geom_boxplot(width=.25, fill='white')
  
  # change axis labels and remove legend:
  p <- p + xlab('') + ylab('hME') + NoLegend()
  
  # save figure
  graphic_filename <- paste0(results_path_ct,"violion_plot_hMEs_",module_name,".pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p)
  dev.off()
}

visualize_module_gene_expression_as_violins <- function(seurat_obj, level, results_path_ct, modules){
  for (module_name in unique(modules$module)) {
    module_df <- subset(modules, module == module_name)
    # Select genes of interest 
    gene.set <- module_df$gene_name
    # Get mean expression of genes of interest per cell
    mean.exp <- colMeans(x = seurat_obj@assays[["RNA"]]@data[gene.set, ], na.rm = TRUE)
    
    # Add mean expression values in 'object@meta.data$gene.set.score'
    if (all(names(x = mean.exp) == rownames(x = seurat_obj@meta.data))) {
      seurat_obj@meta.data[[module_name]] <- mean.exp
    }
    
    # Plot INH-M4 hME using Seurat VlnPlot function
    p <- VlnPlot(
      seurat_obj,
      features = module_name,
      group.by = level,
      pt.size = 0 # don't show actual data points
    )
    
    # add box-and-whisker plots on top:
    p <- p + geom_boxplot(width=.25, fill='white')
    
    # change axis labels and remove legend:
    p <- p + xlab('') + ylab('Average Gene Expression') + NoLegend()
    
    # save figure
    graphic_filename <- paste0(results_path_ct,"violion_plot_average_gene_expression_",module_name,".pdf")
    pdf(file = graphic_filename,   # The directory you want to save the file in
        width = 12, # The width of the plot in inches
        height = 6) # The height of the plot in inches
    print(p)
    dev.off()
  }
}

visualize_number_cells_per_donor_proportion = function(seurat_obj,opt_data,add_str){
  #plot number of cells per sample per celltype --> find settings such that "enough" is covered 
  #(kicking put some cell types is ok, but metacells should be build on a range of control and scz samples)

  df <- data.frame(donor = seurat_obj@meta.data$Donor, celltype=seurat_obj@meta.data$cluster_name_15CTs)
  df <- count(df[,])
  CT_list <- unique(df$celltype)
  target_proportion_samples_vec <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)
  for (ct in CT_list){
    df_sel_CT <- df[df$celltype==ct, ]
    df_sel_CT_sorted <- df_sel_CT[order(df_sel_CT$freq,decreasing=TRUE),]
    total_number_donors_CT <- dim(df_sel_CT)[1]
    df_R <- data.frame(celltype=ct,target_proportion_samples = target_proportion_samples_vec,number_of_cells=0)
    for (target_proportion_samples in target_proportion_samples_vec){
      number_donors_to_include <- floor(target_proportion_samples*total_number_donors_CT)
      i <- which(target_proportion_samples_vec == target_proportion_samples)
      if (number_donors_to_include==0){
        df_R$number_of_cells[i] <- 0
      } else{
        df_R$number_of_cells[i] <- min(df_sel_CT_sorted$freq[1:number_donors_to_include])
      }
    }
    if (ct==CT_list[1]){
      DF_R <- df_R
    } else{
      DF_R <-rbind(DF_R,df_R)
    }
  }
  
  #save  result
  write.csv(DF_R, paste0(opt_data$results_path,"number_cells_per_donor",add_str,".csv"))
  
  #visualize result
  p <- ggplot(data=DF_R, aes(x=target_proportion_samples, y=number_of_cells, group = celltype)) +
    geom_line(aes(color=celltype)) + 
    geom_point(aes(color=celltype,linetype=celltype)) +
    theme_minimal()
  
  graphic_filename <- paste0(opt_data$results_path,"metacell_exploration_min_number_cells",add_str,".pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 9) # The height of the plot in inches
  print(p)
  dev.off()
  
  return(DF_R)
}

visualize_DMEs_as_lollipop = function(seurat_obj, DMEs, settings, ct, ct_, results_path, opt_test){
  # plot DMEs - positive log fold change is over-expressed in group 1
  # the size of each dot corresponds to the number of genes in that module.
  # 'X' over point indicates non significance 
  #browser()
  p1 <- PlotDMEsLollipop(
    seurat_obj, 
    DMEs, 
    wgcna_name=ct, 
    pvalue = "p_val_adj"
  )
  
  graphic_filename_p1 <- paste0(results_path,"DME_lollipop_", ct_,"_", opt_test, ".pdf")
  pdf(file = graphic_filename_p1,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  print(p1)
  dev.off()
}

visualize_DMEs_as_volcano = function(seurat_obj, DMEs, settings, ct, ct_, results_path, opt_test) {
  # plot DMEs - positive log fold change is over-expressed in group 1
  # dots in grey zone are not significant 
  # only significant modules are labelled 
  p2 <- PlotDMEsVolcano(seurat_obj,
                        DMEs,
                        plot_labels=TRUE,
                        label_size=2,
                        wgcna_name = ct) + xlim(min(DMEs$avg_log2FC),max(DMEs$avg_log2FC))
  
  graphic_filename_p2 <- paste0(results_path,"DME_volcano_", ct_, "_", opt_test,".pdf")
  pdf(file = graphic_filename_p2,   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(p2)
  dev.off()
}

visualize_network <- function(seurat_obj, results_path){
  #browser()
  #Module network plot:
  #each node represents a gene, and each edge represents the co-expression relationship between two genes in the network
  #top 10 hub genes by kME are placed in the center of the plot, while the remaining 15 genes are placed in the outer circle
  ModuleNetworkPlot_v2(seurat_obj, 
                          plot_size = c(6, 6),
                          edge.alpha = 0.3, 
                          vertex.size = 2, 
                          vertex.label.cex=1,
                          outdir = paste0(results_path,"/Module_Network_Plot/"))
  
  
  
  # hubgene network (all modules together)
  # takes the top n hub genes as specified by the user, and other randomly selected genes
  # number of edges in the network can be downsampled using the edge_prop parameter
  # each node represents a gene and each edge represents a co-expression relationship. 
  # color intramodular edges: the module's color
  # intermodular edges: gray
  graphic_filename_p1 <- paste0(results_path,"Hub_Gene_network.pdf")
  pdf(file = graphic_filename_p1,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=22,
    edge_prop = 0.5,
    mods = 'all'
  )
  dev.off()
  
  ##visualize several genes simultaneously with UMAP
  #UMAP of topological overlap matrix (TOM)
  #specify number of hub genes per module 
  #unsupervised:
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10,
    n_neighbors = 15, #neighbors parameter for UMAP
    min_dist = 0.1 # min distance between points in UMAP space
  )
  #get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)
  #plot with ggplot
  p2_unsupervised <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(
      color = umap_df$color, # color each point by WGCNA module
      size=umap_df$kME*2 #size of each point based on intramodular connectivity
    ) +
    umap_theme()
  #save plot
  graphic_filename_p2_unsupervised <- paste0(results_path,"Hub_Gene_network_UMAP_unsupervised.pdf")
  pdf(file = graphic_filename_p2_unsupervised,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 9) # The height of the plot in inches
  print(p2_unsupervised)
  dev.off()
  
  ##visualize genes and their co-expression relationships
  #save
  graphic_filename_p3 <- paste0(results_path,"Co_expression_network_UMAP.pdf")
  pdf(file = graphic_filename_p3,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.3,
    sample_edges = TRUE,
    edge_prop = 0.1,#proportion of edges to sample (10%)
    label_hubs = 3, # how many hub genes should be labeled per module
    keep_grey_edges = FALSE,
    vertex.label.cex=0.1
  )
  dev.off()
  
  ##visualize several genes simultaneously with UMAP
  #UMAP of topological overlap matrix (TOM)
  #specify number of hub genes per module 
  #supervised: using the infor to which module a gene belongs
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10,
    n_neighbors = 15, #neighbors parameter for UMAP
    min_dist = 0.1, # min distance between points in UMAP space
    supervised=TRUE,
    target_weight = 0.5 # weight close to 0 almost unsupervised (preserves structure of the data), weight close to 1 weights heavily based on the module labels
  )
  #get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)
  #plot with ggplot
  p2_supervised <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(
      color = umap_df$color, # color each point by WGCNA module
      size=umap_df$kME*2 #size of each point based on intramodular connectivity
    ) +
    umap_theme()
  #save plot
  graphic_filename_p2_supervised <- paste0(results_path,"Hub_Gene_network_UMAP_supervised.pdf")
  pdf(file = graphic_filename_p2_supervised,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 9) # The height of the plot in inches
  print(p2_supervised)
  dev.off()
  
  return(seurat_obj)
}

visualize_significant_modules_as_networks <- function(seurat_obj, DMEs, results_path){
  sign_modules <- deframe(DMEs[DMEs$p_val_adj<=0.05,"module"])
  for (module_name in sign_modules){
    ## hub gene network
    graphic_filename <- paste0(results_path,"Module_Network_Plot/","Hub_Gene_network_module",module_name,".pdf")
    pdf(file = graphic_filename,   # The directory you want to save the file in
        width = 6, # The width of the plot in inches
        height = 6) # The height of the plot in inches
    HubGeneNetworkPlot(
      seurat_obj,
      n_hubs = 15, 
      n_other=22,#change to all
      edge_prop = 0.5,
      mods = module_name #change to module_name
    )
    dev.off()
    
    ## functional networks

    
  }
}


#add trycatch:
ModuleNetworkPlot_v2 <- function (seurat_obj, mods = "all", outdir = "ModuleNetworks", 
          plot_size = c(6, 6), wgcna_name = NULL, label_center = FALSE, 
          edge.alpha = 0.25, vertex.label.cex = 1, vertex.size = 6, 
          ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  if (mods == "all") {
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
  }
  if (!all(paste0("kME_", as.character(mods)) %in% colnames(modules))) {
    stop("Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?")
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  cat(paste0("Writing output files to ", outdir))
  TOM <- GetTOM(seurat_obj, wgcna_name)
  n_hubs <- 25
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c("gene_name", paste0("kME_", cur_mod))] %>% 
      top_n(n_hubs)
    colnames(cur)[2] <- "var"
    cur %>% arrange(desc(var)) %>% .$gene_name
  })
  names(hub_list) <- mods
  for (cur_mod in mods) {
    print(cur_mod)
    c <-tryCatch(
      { plot_network(modules,
                     cur_mod,
                     TOM,
                     label_center,
                     outdir,
                     plot_size,
                     edge.alpha,
                     vertex.label.cex,
                     vertex.size,
                     ...)
      },
      error=function(cond){
        message("ModuleNetworkPlot failed!")
        message(cond)
        return(c())
      }
    )
  }
}

plot_network <- function(modules,cur_mod,TOM,label_center,outdir,plot_size,edge.alpha,vertex.label.cex,vertex.size,...){

  cur_color <- modules %>% subset(module == cur_mod) %>% 
    .$color %>% unique
  n_genes = 25
  n_conns = 500
  cur_kME <- paste0("kME_", cur_mod)
  cur_genes <- hub_list[[cur_mod]]
  matchind <- match(cur_genes, colnames(TOM))
  reducedTOM = TOM[matchind, matchind]
  orderind <- order(reducedTOM, decreasing = TRUE)
  connections2keep <- orderind[1:n_conns]
  reducedTOM <- matrix(0, nrow(reducedTOM), ncol(reducedTOM))
  reducedTOM[connections2keep] <- 1
  if (label_center) {
    cur_genes[11:25] <- ""
  }
  gA <- graph.adjacency(as.matrix(reducedTOM[1:10, 1:10]), 
                        mode = "undirected", weighted = TRUE, diag = FALSE)
  gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes, 
                                             11:n_genes]), mode = "undirected", weighted = TRUE, 
                        diag = FALSE)
  layoutCircle <- rbind(layout.circle(gA)/2, layout.circle(gB))
  g1 <- graph.adjacency(as.matrix(reducedTOM), mode = "undirected", 
                        weighted = TRUE, diag = FALSE)
  pdf(paste0(outdir, "/", cur_mod, ".pdf"), width = plot_size[1], 
      height = plot_size[2], useDingbats = FALSE)
  plot(g1, edge.color = adjustcolor(cur_color, alpha.f = 0.25), 
       edge.alpha = edge.alpha, vertex.color = cur_color, 
       vertex.label = as.character(cur_genes), vertex.label.dist = 1.1, 
       vertex.label.degree = -pi/4, vertex.label.color = "black", 
       vertex.label.family = "Helvetica", vertex.label.font = 3, 
       vertex.label.cex = vertex.label.cex, vertex.frame.color = "black", 
       layout = jitter(layoutCircle), vertex.size = vertex.size, 
       main = paste(cur_mod))
  dev.off()  
}

visualize_module_trait_corr <- function(seurat_obj, results_path){
  # heatmap of correlation between traits and modules
  p <- PlotModuleTraitCorrelation_v2(
    seurat_obj,
    label = 'fdr',
    label_symbol = 'stars',
    text_size = 1,
    text_digits = 2,
    text_color = 'black',
    high_color = 'darkcyan',
    mid_color = 'white',
    low_color = 'darkslateblue',
    plot_max = 0.2,
    combine=TRUE
  )
  graphic_filename <- paste0(results_path,"module_trait_correlation_heatmap.pdf")
  pdf(file = graphic_filename,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  print(p)
  dev.off()
}

PlotModuleTraitCorrelation_v2 <- function (seurat_obj, high_color = "red", mid_color = "grey90", 
                                           low_color = "blue", label = NULL, label_symbol = "stars", 
                                           plot_max = NULL, text_size = 2, text_color = "black", text_digits = 3, 
                                           combine = TRUE, wgcna_name = NULL){
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  temp <- GetModuleTraitCorrelation(seurat_obj, wgcna_name)
  cor_list <- temp$cor
  pval_list <- temp$pval
  fdr_list <- temp$fdr
  if (is.null(dim(cor_list[[1]]))) {
    stop("ModuleTraitCorrelation was run only for one trait. Heatmaps are not suggested for visualizing only one variable!")
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  module_colors <- modules %>% dplyr::select(c(module, color)) %>% 
    distinct %>% subset(module != "grey") %>% arrange(module)
  mod_colors <- module_colors$color
  module_colors$var <- 1
  module_colorbar <- module_colors %>% ggplot(aes(x = module, 
                                                  y = var, fill = module)) + geom_tile() + scale_fill_manual(values = mod_colors) + 
    NoLegend() + RotatedAxis() + theme(plot.title = element_blank(), 
                                       axis.line = element_blank(), axis.ticks.y = element_blank(), 
                                       axis.text.y = element_blank(), axis.title = element_blank(), 
                                       plot.margin = margin(0, 0, 0, 0))
  plot_list <- list()
  for (i in names(cor_list)) {
    cor_mat <- as.matrix(cor_list[[i]])
    pval_mat <- as.matrix(pval_list[[i]])
    fdr_mat <- as.matrix(fdr_list[[i]])
    print(i)
    plot_df <- reshape2::melt(cor_mat)
    colnames(plot_df) <- c("Trait", "Module", "cor")
    if (!is.null(label)) {
      if (label == "fdr") {
        p_df <- reshape2::melt(fdr_mat)
      }
      else if (label == "pval") {
        p_df <- reshape2::melt(pval_mat)
      }
      #browser()
      #last column should be label
      colnames(p_df) <- c("Trait", "Module", "pval")
      plot_df$pval <- p_df$pval
      print(levels(plot_df$Trait))
      #browser()
      if (label_symbol == "stars") {
        plot_df$significance <- gtools::stars.pval(plot_df$pval)
      }
      else if (label_symbol == "numeric") {
        plot_df$significance <- ifelse(plot_df$pval <= 
                                         0.05, formatC(plot_df$pval, digits = text_digits), 
                                       "")
      }
      else {
        stop("Invalid input for label_symbol. Valid choices are stars or numeric.")
      }
    }
    if (is.null(plot_max)) {
      max_plot <- max(abs(range(plot_df$cor)))
    }
    else {
      max_plot <- plot_max
      plot_df$cor <- ifelse(abs(plot_df$cor) >= plot_max, 
                            plot_max * sign(plot_df$cor), plot_df$cor)
    }
    p <- ggplot(plot_df, aes(x = Module, y = as.numeric(Trait), 
                             fill = cor)) + geom_tile() + scale_fill_gradient2(limits = c(-1 * 
                                                                                            max_plot, max_plot), high = high_color, mid = mid_color, 
                                                                               low = low_color, guide = guide_colorbar(ticks = FALSE, 
                                                                                                                       barwidth = 16, barheight = 0.5)) + scale_y_continuous(breaks = 1:length(levels(plot_df$Trait)), 
                                                                                                                                                                             labels = levels(plot_df$Trait), sec.axis = sec_axis(~., 
                                                                                                                                                                                                                                 breaks = 1:length(levels(plot_df$Trait)), labels = levels(plot_df$Trait))) + 
      RotatedAxis() + ylab("") + xlab("") + ggtitle(i) + 
      theme(plot.title = element_text(hjust = 0.5), axis.line = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y.left = element_blank(), 
            axis.title.y = element_text(angle = 0, vjust = 0.5), 
            legend.title = element_blank(), legend.position = "bottom")
    if (!is.null(label)) {
      p <- p + geom_text(label = plot_df$significance, 
                         color = text_color, size = text_size)
    }
    plot_list[[i]] <- p
  }
  if (combine) {
    for (i in 1:length(plot_list)) {
      plot_list[[i]] <- plot_list[[i]] + ylab(names(plot_list)[i]) + 
        theme(plot.margin = margin(t = 0, r = 0, b = 0, 
                                   l = 0), axis.title.x = element_blank(), plot.title = element_blank(), 
              legend.position = "bottom", axis.text.x = element_blank(), 
              axis.ticks = element_blank(), axis.title.y = element_text(angle = 0, 
                                                                        vjust = 0.5))
    }
    plot_list[["module"]] <- module_colorbar
    out <- wrap_plots(plot_list, ncol = 1) + plot_layout(guides = "collect", 
                                                         heights = c(rep(1, length(plot_list) - 1), 0.15)) + 
            plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5), 
                                          legend.position = "bottom", legend.justification = 0.5))
    return(out)
  }
  else {
    return(plot_list)
  }
}

visualize_marker_overlap <- function(overlap_df, results_path, ct){
  # overlap barplot, produces a plot for each cell type
  plot_list <- OverlapBarPlot_v2(overlap_df, fontsize_title = 10)
  graphic_filename_p1 <- paste0(results_path,"overlap_modules_marker_genes_barplot.pdf")
  pdf(file = graphic_filename_p1,   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 14) # The height of the plot in inches
  # stitch plots with patchwork
  print(wrap_plots(plot_list, ncol=3))
  dev.off()
  
  # plot odds ratio of the overlap as a dot plot
  p <- OverlapDotPlot(
    overlap_df,
    plot_var = 'odds_ratio') +
    ggtitle(paste0('Overlap of ', ct, ' modules & cell type markers'))
  
  graphic_filename_p2 <- paste0(results_path,"overlap_modules_marker_genes_dotplot.pdf")
  pdf(file = graphic_filename_p2,   # The directory you want to save the file in
      width = 9, # The width of the plot in inches
      height = 12) # The height of the plot in inches
  print(p)
  dev.off()
}

OverlapBarPlot_v2 <- function (overlap_df, fontsize_title = 10, plot_var = "odds_ratio", logscale = FALSE, neglog = FALSE, label_size = 2, ...) {
  label <- plot_var
  if (plot_var == "odds_ratio") {
    yint <- 1
  }
  else if (plot_var == "fdr") {
    yint <- 0.05
  }
  if (logscale) {
    overlap_df[[plot_var]] <- log(overlap_df[[plot_var]])
    label <- paste0("log(", plot_var, ")")
    yint = log(yint)
  }
  if (neglog) {
    overlap_df[[plot_var]] <- -1 * log(overlap_df[[plot_var]])
    label <- paste0("-log(", label, ")")
    yint = -1 * log(yint)
  }
  groups <- overlap_df$group %>% as.character %>% unique
  plot_list <- list()
  for (cur_group in groups) {
    cur_df <- overlap_df %>% subset(group == cur_group)
    p <- cur_df %>% ggplot(aes(x = reorder(module, get(plot_var)), 
                               y = get(plot_var))) + geom_bar(stat = "identity", 
                                                              fill = cur_df$color) + coord_flip() + xlab("") + 
      ylab(label) + ggtitle(cur_group) + theme(axis.line.y = element_blank(), 
                                               axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                                               plot.title = element_text(hjust = 0.5, size=fontsize_title))
    if (plot_var == "fdr" | plot_var == "odds_ratio") {
      p <- p + geom_hline(yintercept = yint, linetype = "dashed", 
                          color = "gray")
    }
    p <- p + geom_text(aes(label = module, x = module, y = get(plot_var)), 
                       color = "black", size = label_size, hjust = "inward")
    plot_list[[cur_group]] <- p
  }
  plot_list
}


ModuleCorrelogram_v2 <- function (seurat_obj, MEs2 = NULL, features = "hMEs", order = "original", 
                                  method = "ellipse", exclude_grey = TRUE, type = "upper", 
                                  tl.col = "black", tl.srt = 45, sig.level = c(1e-04, 0.001, 
                                                                               0.01, 0.05), pch.cex = 0.7, col = NULL, ncolors = 200, 
                                  wgcna_name = NULL, wgcna_name2 = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (features == "hMEs") {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  }
  else if (features == "MEs") {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  }
  else if (features == "scores") {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  }
  else if (features == "average") {
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  }
  else (stop("Invalid feature selection. Valid choices: hMEs, MEs, scores, average"))
  MEs <- as.matrix(MEs)
  if (exclude_grey) {
    MEs <- MEs[, colnames(MEs) != "grey"]
  }
  if (is.null(col)) {
    colfunc <- grDevices::colorRampPalette(c("seagreen", 
                                             "white", "darkorchid1"))
    col = colfunc(ncolors)
  }
  if (is.null(MEs2)) {
    res <- Hmisc::rcorr(x = MEs)
  }
  else {
    d1_names <- colnames(MEs)
    d2_names <- colnames(MEs2)
    colnames(MEs) <- paste0(d1_names, "_D1")
    colnames(MEs2) <- paste0(d2_names, "_D2")
    res <- Hmisc::rcorr(x = MEs, y = as.matrix(MEs2))
    res$r <- res$r[!grepl("_D1", colnames(res$r)), grepl("_D1", 
                                                         colnames(res$r))]
    colnames(res$r) <- d1_names
    rownames(res$r) <- d2_names
    res$P <- res$P[!grepl("_D1", colnames(res$P)), grepl("_D1", 
                                                         colnames(res$P))]
    colnames(res$P) <- d1_names
    rownames(res$P) <- d2_names
  }
  res$P[is.na(res$P)] <- 0
  corrplot::corrplot(res$r, p.mat = res$P, type = type, order = order, 
                     method = method, tl.col = tl.col, tl.srt = tl.srt, sig.level = sig.level, 
                     pch.cex = pch.cex, col = col, pch = "n", ...)
}

visualize_cell_type_UMAP <- function(seurat_obj, level, results_path){
  p <- DimPlot(seurat_obj, group.by = level)
  pdf(paste0(results_path, "reference_UMAP_annotated_cell_types.pdf"),
      width = 10)
  print(p)
  dev.off()
} 

plot_dendro_heatmap <- function(seurat_obj, results_path, ct){
  # calculate dissimilarity TOM
  TOM <- GetTOM(seurat_obj)
  dissTOM = 1-TOM
  
  # transform dissTOM to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM^7
  # set the diagonal to NA to improve the clarity of the plot
  diag(plotTOM) = NA
  
  geneTree = seurat_obj@misc[[ct]][["wgcna_net"]][["dendrograms"]][[1]]
  moduleColors <- seurat_obj@misc[[ct]][["wgcna_net"]][["colors"]]
  # plot heatmap
  pdf(paste0(results_path, 'network_heatmap_plot_all_genes.pdf'), 
      width = 10, 
      height = 10)
  TOMplot(plotTOM, geneTree, moduleColors, main = paste0(ct, " network heatmap plot"))
  dev.off()
}