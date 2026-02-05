#about: utils for cell-cell communication inference with cellChat v2
#author: Lisa Bast
#date: 11.03.25
#version: 0.0.1

load_data_and_get_cell_chat_object <- function(data_path,data_file,subset_condition){
  #load a scRNA-seq data matrix (gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames)
  #and its associated cell meta data
  data <- loomR::connect(paste0(data_path,data_file),mode = "r", skip.validate = TRUE)
  
  gene_IDs <- data[['row_attrs/Accession']][]
  genes <- data[['row_attrs/Gene']][]
  cell_groups <- data[['col_attrs/cluster_name_15CTs']][]
  
  #data matrix selection without removed nuclei:
  idx <- which(!(cell_groups %in% c('none (removed)')))
  cellIDs <- data[['col_attrs/CellID']][idx]
  condition <- data[['col_attrs/Disease']][idx]
  cell_type <- cell_groups[idx]
  cell_type <- simplify_cell_type_names(cell_type)
  donor <- data[['col_attrs/Donor']][idx]
  counts <- t(data[["matrix"]][idx,])
  rownames(counts) <- genes
  colnames(counts) <- cellIDs
  #normalize data (library-size normalization and then log-transformed with a pseudocount of 1)
  counts_normalized <- normalizeData(counts)
  
  data$close_all()
  
  # a dataframe with rownames containing cell meta data
  meta <- data.frame(condition =  condition,labels = cell_type,donor = donor)
  row.names(meta) <- cellIDs
  if (subset_condition=="both"){
    cell.use = rownames(meta)
  } else{
    cell.use = rownames(meta)[meta$condition == subset_condition] # extract the cell names from disease data
  }
  
  # Prepare input data for CelChat analysis
  counts_normalized = counts_normalized[, cell.use]
  meta = meta[cell.use, ]
  # meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
  unique(meta$labels) # check the cell labels
  
  cellchat <- createCellChat(object = counts_normalized, meta = meta, group.by = "labels")
  
  return(cellchat)
}

simplify_cell_type_names <- function(cell_type){
  #replace Inhibitory with Inh, Excitatory with Exc, Layer with L, neurons with " "
  cell_type_simplified <- str_replace_all(cell_type,"Excitatory","Exc")
  cell_type_simplified <- str_replace_all(cell_type_simplified,"Inhibitory","Inh")
  cell_type_simplified <- str_replace_all(cell_type_simplified," neurons","")
  cell_type_simplified <- str_replace_all(cell_type_simplified,"Layer ","L")
  cell_type_simplified <- str_replace_all(cell_type_simplified,"Microglial cells","Microglia")
  cell_type_simplified <- str_replace_all(cell_type_simplified,"Endothelial and mural cells","Endothelial/mural")
  cell_type_simplified <- str_replace_all(cell_type_simplified,"Oligodendrocyte progenitor cells","OPCs")
  return(cell_type_simplified)
}

visualize_aggr_network <- function(cellchat_combined,cellchat_SCZ,cellchat_CTRL,results_path){
  ## visuzalize the aggregated network
  
  #For example, showing the number of interactions or the total interaction strength (weights) 
  #between any two cell groups using circle plot.
  pdf(file = paste0(results_path,"aggregated_network_number_interactions_and_interaction_strength.pdf"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  colors <- get_ct_colors(cellchat_combined)
  groupSize_SCZ <- as.numeric(table(cellchat_SCZ@idents))
  groupSize_CTRL <- as.numeric(table(cellchat_CTRL@idents))
  par(mfrow = c(2,2), xpd=TRUE)
  netVisual_circle(cellchat_SCZ@net$count, vertex.weight = groupSize_SCZ, vertex.label.cex=0.4, weight.scale = T, label.edge= F, title.name = "Number of interactions SCZ",color.use = colors)
  netVisual_circle(cellchat_CTRL@net$count, vertex.weight = groupSize_CTRL, vertex.label.cex=0.4, weight.scale = T, label.edge= F, title.name = "Number of interactions CTRL",color.use = colors)
  netVisual_circle(cellchat_SCZ@net$weight, vertex.weight = groupSize_SCZ, vertex.label.cex=0.4, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength SCZ",color.use = colors)
  netVisual_circle(cellchat_CTRL@net$weight, vertex.weight = groupSize_CTRL, vertex.label.cex=0.4, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength CTRL",color.use = colors)
  dev.off()
  
  #Due to the complicated cell-cell communication network, we can examine the signaling sent from 
  #each cell group. Here we also control the parameter edge.weight.max so that we can compare edge 
  #weights between differet networks. 
  pdf(file = paste0(results_path,"aggregated_network_signaling_sent_from_each_celltype_SCZ.pdf"),   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 15) # The height of the plot in inches
  
  mat_SCZ <- cellchat_SCZ@net$weight
  par(mfrow = c(4,4), xpd=TRUE)
  for (i in 1:nrow(mat_SCZ)) {
    mat2 <- matrix(0, nrow = nrow(mat_SCZ), ncol = ncol(mat_SCZ), dimnames = dimnames(mat_SCZ))
    mat2[i, ] <- mat_SCZ[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize_SCZ, vertex.label.cex = 0.7, weight.scale = T, edge.weight.max = max(mat_SCZ), title.name = rownames(mat_SCZ)[i])
  }
  dev.off()
  
  pdf(file = paste0(results_path,"aggregated_network_signaling_sent_from_each_celltype_CTRL.pdf"),   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 15) # The height of the plot in inches
  
  mat_CTRL <- cellchat_CTRL@net$weight
  par(mfrow = c(4,4), xpd=TRUE)
  for (i in 1:nrow(mat_CTRL)) {
    mat2 <- matrix(0, nrow = nrow(mat_CTRL), ncol = ncol(mat_CTRL), dimnames = dimnames(mat_CTRL))
    mat2[i, ] <- mat_CTRL[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize_CTRL, vertex.label.cex = 0.7, weight.scale = T, edge.weight.max = max(mat_CTRL), title.name = rownames(mat_CTRL)[i])
  }
  dev.off()
}

visualize_cell_cell_comm_network <- function(cellchat,results_path,pathways.show,group_str){
  ##Visualize the cell-cell communication network
  #several ways for visualizing cell-cell communication network, 
  #including hierarchical plot, circle plot, Chord diagram, and bubble plot.
  
  #a) Chord diagram
  pdf(file = paste0(results_path,"Chord_diagram_",group_str,".pdf"),   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 15) # The height of the plot in inches
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  
  
  # b) Heatmap of communication probabilities
  pdf(file = paste0(results_path,"Heatmap_communication_probabilities_",group_str,".pdf"),   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat,signaling = pathways.show, color.heatmap = "Reds"))
  dev.off()

}

visualize_number_of_interactions_and_interaction_strength <- function(cellchat_combined, results_path){
  gg1 <- compareInteractions(cellchat_combined, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat_combined, show.legend = F, group = c(1,2), measure = "weight")
  gg<-gg1 + gg2
  ggsave(paste0(results_path,"number_of_interactions_and_interaction_strength_comparison.pdf"),gg,width = 5, height = 5)
}

visualize_differential_number_of_interactions_and_interaction_strength <- function(cellchat_combined, results_path,opt_heatmap){
  #"count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
  #turquoise (blue): decreased signalling in SCZ
  #darkred (red): increased signalling in SCZ
  pdf(file = paste0(results_path,"differential_number_of_interactions_and_interaction_strength_comparison_circle_plot.pdf"),   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 15) # The height of the plot in inches
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat_combined,measure="count",comparison = c(1,2),color.edge = c("turquoise4","darkred"),weight.scale = T)
  netVisual_diffInteraction(cellchat_combined,measure="weight",comparison = c(1,2),color.edge = c("turquoise4","darkred"),weight.scale = T)
  dev.off()
  
  if (opt_heatmap == "pathways"){
    slot_name = "netP"
  }else{ #ligand-receptor pair
    slot_name = "net"
  }
  
  pdf(file=paste0(results_path,"heatmap_differential_number_of_interactions_",opt_heatmap,".pdf"),
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(1,2))
  print(netVisual_heatmap(cellchat_combined, slot.name = slot_name, color.heatmap = c("darkred","turquoise4"),font.size = 10,font.size.title = 14))
  dev.off()
  
  pdf(file=paste0(results_path,"heatmap_differential_interaction_strength_",opt_heatmap,".pdf"),
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(1,1))
  print(netVisual_heatmap(cellchat_combined, slot.name = slot_name, color.heatmap = c("darkred","turquoise4"), measure = "weight",font.size = 10,font.size.title = 14))
  dev.off()
}

visualize_ingoing_vs_outgoing_interaction_strength_per_celltype <- function(cellchat_SCZ,cellchat_CTRL,results_path){
  #Compare the major sources and targets in 2D space
  # identification of the cell populations with significant changes in sending or receiving signals between different datasets. 
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  
  colors <- get_ct_colors(cellchat_combined,with_white=FALSE)

  gg <- list()
  gg[[1]] <- netAnalysis_signalingRole_scatter(cellchat_SCZ, title = "SCZ", weight.MinMax = weight.MinMax,label.size = 5,font.size = 16,font.size.title = 14, dot.alpha = 0.5,color.use = colors)
  gg[[2]] <- netAnalysis_signalingRole_scatter(cellchat_CTRL, title = "CTRL", weight.MinMax = weight.MinMax,label.size = 5,font.size = 16,font.size.title = 14,color.use = colors)
  pdf(file=paste0(results_path,"ingoing_vs_outgoing_interaction_strength_per_celltype.pdf"),
      width = 15, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(1,1))
  print(patchwork::wrap_plots(plots = gg, axes="collect_y"))
  dev.off()
  
  colors <- get_ct_colors(cellchat_combined,with_white=TRUE)
  gg <- netAnalysis_signalingRole_scatter_comparison(cellchat_CTRL,cellchat_SCZ, title = "CTRL vs. SCZ", weight.MinMax = weight.MinMax,label.size = 5,font.size = 16,font.size.title = 14, color.use = colors) 

  pdf(file=paste0(results_path,"ingoing_vs_outgoing_interaction_strength_per_celltype_combined.pdf"),
      width = 15, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(1,1))
  print(gg)
  dev.off()
}

netAnalysis_signalingRole_scatter_comparison = function (object1,object2, signaling = NULL, color.use = NULL, slot.name = "netP", 
          group = NULL, weight.MinMax = NULL, dot.size = c(2, 6), 
          point.shape = c(21, 22, 24, 23, 25, 8, 3), label.size = 3, 
          dot.alpha = 0.6, x.measure = "outdeg", y.measure = "indeg", 
          xlabel = "Outgoing interaction strength", ylabel = "Incoming interaction strength", 
          title = NULL, font.size = 10, font.size.title = 10, do.label = T, 
          show.legend = T, show.axes = T) 
{
  if ((length(slot(object1, slot.name)$centr) == 0)|(length(slot(object2, slot.name)$centr) == 0)) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  if (sum(c(x.measure, y.measure) %in% names(slot(object1, 
                                                  slot.name)$centr[[1]])) != 2) {
    stop(paste0("`x.measure, y.measure` should be one of ", 
                paste(names(slot(object, slot.name)$centr[[1]]), 
                      collapse = ", "), "\n", "`outdeg_unweighted` is only supported for version >= 1.1.2"))
  }
  if (sum(c(x.measure, y.measure) %in% names(slot(object2, 
                                                  slot.name)$centr[[1]])) != 2) {
    stop(paste0("`x.measure, y.measure` should be one of ", 
                paste(names(slot(object, slot.name)$centr[[1]]), 
                      collapse = ", "), "\n", "`outdeg_unweighted` is only supported for version >= 1.1.2"))
  }

  df1 <- get_df(object1,slot.name,signaling,group,x.measure,y.measure)
  df2 <- get_df(object2,slot.name,signaling,group,x.measure,y.measure)
  
  df1$fill_color <- df1$labels
  df2$fill_color <- "white"
  df <- rbind(df1,df2)

  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(object1@idents)+1)
  }
  if (!is.null(group)) {
    gg <- ggplot(data = df, aes(x, y)) + 
          geom_point(aes(size = Count, 
                         colour = labels, 
                         fill = fill_color, 
                         shape = Group))
  }
  else {
    gg <- ggplot(data = df, aes(x, y)) + 
          geom_point(aes(size = Count, 
                         colour = labels, 
                         fill = fill_color))
  }
  gg <- gg + CellChat_theme_opts() + theme(text = element_text(size = font.size), 
                                           legend.key.height = grid::unit(0.15, "in")) + labs(title = title, 
                                                                                              x = xlabel, y = ylabel) + theme(plot.title = element_text(size = font.size.title, 
                                                                                                                                                        face = "plain")) + theme(axis.line.x = element_line(size = 0.25), 
                                                                                                                                                                                 axis.line.y = element_line(size = 0.25))
  #gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
  #                                                     alpha = dot.alpha), drop = FALSE) + guides(fill = FALSE)
  gg <- gg + scale_fill_manual(values = factor(color.use)) 
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE) + 
                                guides(colour = FALSE)
  if (!is.null(group)) {
    gg <- gg + scale_shape_manual(values = point.shape[1:length(unique(df$Group))])
  }
  if (is.null(weight.MinMax)) {
    gg <- gg + scale_size_continuous(range = dot.size)
  }
  else {
    gg <- gg + scale_size_continuous(limits = weight.MinMax, 
                                     range = dot.size)
  }
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(mapping = aes(label = labels, 
                                                      colour = labels), size = label.size, show.legend = F, 
                                        segment.size = 0.2, segment.alpha = 0.5)
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  gg
}

get_df <- function(object,slot.name,signaling,group,x.measure,y.measure){
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]][[x.measure]]
    incoming[, i] <- centr[[i]][[y.measure]]
  }
  if (is.null(signaling)) {
    message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
  }
  else {
    message("Signaling role analysis on the cell-cell communication network from user's input")
    signaling <- signaling[signaling %in% object@netP$pathways]
    if (length(signaling) == 0) {
      stop("There is no significant communication for the input signaling. All the significant signaling are shown in `object@netP$pathways`")
    }
    outgoing <- outgoing[, signaling, drop = FALSE]
    incoming <- incoming[, signaling, drop = FALSE]
  }
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  num.link <- aggregateNet(object, signaling = signaling, 
                           return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link) - diag(num.link)
  df <- data.frame(x = outgoing.cells, y = incoming.cells, 
                   labels = names(incoming.cells), Count = num.link)
  df$labels <- factor(df$labels, levels = names(incoming.cells))
  if (!is.null(group)) {
    df$Group <- group
  }  
  return(df)
}

visualize_ingoing_vs_outgoing_signaling_changes_per_celltype <- function(cellchat_combined,results_path){
  CTs <- levels(cellchat_combined@idents$joint)
  for (ct in CTs){
    ct_str <- adapt_celltype_str(ct)
    pdf(file=paste0(results_path,"ingoing_vs_outgoing_interaction_signaling_change_",ct_str,".pdf"),
        width = 10, # The width of the plot in inches
        height = 10) # The height of the plot in inches
    par(mfrow=c(1,1))
    print(netAnalysis_signalingChanges_scatter(cellchat_combined, idents.use = ct, do.label=T, font.size = 12, font.size.title = 15))
    dev.off()
  }
}

adapt_celltype_str <- function(ct){
  ct_str <- gsub("/","_",ct)
  ct_str <- gsub(" ","_",ct_str)
  return(ct_str)
}

visualize_ranked_signaling_pathways <- function(cellchat_combined, measure, source_celltype, target_celltype, results_path){
  title_str <- paste0("signaling pathways from ",source_celltype," to ",target_celltype)
  source_celltype_str <- adapt_celltype_str(source_celltype)
  target_celltype_str <-  adapt_celltype_str(target_celltype)
  gg1 <- rankNet(cellchat_combined, measure=measure, mode = "comparison", sources.use = source_celltype, targets.use = target_celltype, stacked = T, do.stat = TRUE, color.use=c("steelblue","darkorange"))
  gg2 <- rankNet(cellchat_combined, measure =measure, mode = "comparison", sources.use = source_celltype, targets.use = target_celltype, stacked = F, do.stat = TRUE, color.use=c("steelblue","darkorange"))
  gg <- gg1 + gg2
  
  pdf(file=paste0(results_path,"ranked_signaling_pathways_from_",source_celltype_str,"_to_",target_celltype_str,"_ranked_by_",measure,".pdf"),
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(3,1))
  print(gg + plot_annotation(title = title_str, theme = theme(plot.title = element_text(hjust = 0.5))))
  dev.off()
}

visualize_difference_in_communication_probabilities <- function(cellchat_combined,source_celltype,target_celltype,results_path){
  source_celltype_str <- adapt_celltype_str(source_celltype)
  target_celltype_str <-  adapt_celltype_str(target_celltype)
  #dysfunctional signaling: comparison of communication probabilities
  gg1 <- netVisual_bubble(cellchat_combined, sources.use = source_celltype, targets.use = target_celltype,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in SCZ", angle.x = 45, remove.isolate = T)
  gg2 <- netVisual_bubble(cellchat_combined, sources.use = source_celltype, targets.use = target_celltype,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in SCZ", angle.x = 45, remove.isolate = T)
  gg <- gg1 + gg2
  
  pdf(file=paste0(results_path,"difference_in_communication_probabilities_from_",source_celltype_str,"_to_",target_celltype_str,".pdf"),
      width = 10, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow=c(3,1))
  print(gg)
  dev.off()
}
  
visualize_DEGS_bubble <- function(cellchat_combined,object.list,source_celltype,target_celltype,results_path){

    filename <- paste0(results_path,"DEGs_from_",source_celltype,"_to_",target_celltype,"_bubble.pdf")
      
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    gg1 <- netVisual_bubble(cellchat_combined, pairLR.use = pairLR.use.up, sources.use = source_celltype, targets.use = target_celltype, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    pairLR.use.down = net.down[, "interaction_name", drop = F]
    gg2 <- netVisual_bubble(cellchat_combined, pairLR.use = pairLR.use.down, sources.use = source_celltype, targets.use = target_celltype, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    gg <- gg1 + gg2
    pdf(file=filename,
        width = 20, # The width of the plot in inches
        height = 10) # The height of the plot in inches
    par(mfrow=c(1,1))
    print(gg)
    dev.off()
}


visualize_DEGS_chord_cell <- function(cellchat_combined,object.list,results_path,net.up,net.down,pw,str_up_down){
    
    #CTs <- levels(cellchat_combined@idents$joint)
    filename <- paste0(results_path,"DEGs_chord_cell_",pw,"_",str_up_down,".pdf")
    colors <- get_ct_colors(cellchat_combined)
    # Chord diagram
    pdf(file=filename,
        width = 20, # The width of the plot in inches
        height = 10) # The height of the plot in inches
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_chord_cell(object.list[[2]], signaling = pw, title.name = paste0(pw,"_",str_up_down," Signaling network - ", names(object.list)[1]),color.use = colors)
    netVisual_chord_cell(object.list[[1]], signaling = pw, title.name = paste0(pw,"_",str_up_down," Signaling network - ", names(object.list)[2]),color.use = colors)
    dev.off()
}

visualize_DEGS_chord_gene <- function(cellchat_combined,object.list,results_path,net,net_L_or_R,add_str_ligand_or_receptor,opt_filter_top5){

  #filter pairs  
  if (opt_filter_top20){
    if(add_str_ligand_or_receptor=="ligand"){
      top_interactions_up <- get_top5_interactions(net_L_or_R, TRUE,TRUE)
      top_interactions_down <- get_top5_interactions(net_L_or_R, TRUE,FALSE)
      #filter list with pairLR.use
      net.up <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ", pairLR.use = top_interactions_up, ligand.logFC = 0.15)
      net.down <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ", pairLR.use = top_interactions_down, ligand.logFC = -0.25)
    } else{
      top_interactions_up <- get_top5_interactions(net_L_or_R, FALSE, TRUE)
      top_interactions_down <- get_top5_interactions(net_L_or_R, FALSE, FALSE)
      #filter list with pairLR.use
      net.up <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ", pairLR.use = top_interactions_up,receptor.logFC = 0.15)
      net.down <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ", pairLR.use = top_interactions_down,receptor.logFC = -0.15)
    }
  } else{
    if(add_str_ligand_or_receptor=="ligand"){
      # extract the ligand-receptor pairs with upregulated ligands in SCZ
      net.up <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ",ligand.pvalues = 0.05,ligand.logFC = 0.18,ligand.pct.1=0.3)#,thresh = 0.05)
      # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in CTRL, i.e.,downregulated in LS
      net.down <- subsetCommunication(cellchat_combined, net = net, datasets = "CTRL",ligand.pvalues = 0.05, ligand.logFC = 0.25,ligand.pct.1=0.3)#,thresh = 0.05)
    } else{
      # extract the ligand-receptor pairs with upregulated ligands in SCZ
      net.up <- subsetCommunication(cellchat_combined, net = net, datasets = "SCZ",receptor.pvalues = 0.05,receptor.logFC = 0.23,receptor.pct.1=0.3)#,thresh = 0.05)
      # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in CTRL, i.e.,downregulated in LS
      net.down <- subsetCommunication(cellchat_combined, net = net, datasets = "CTRL",receptor.pvalues = 0.05,receptor.logFC = 0.15,receptor.pct.1=0.3)#,thresh = 0.05)
    }
  }
  
  #plot
  filename <- paste0(results_path,"DEGs_chord_gene_",add_str_ligand_or_receptor,"_deregulated.pdf")
  colors <- get_ct_colors(cellchat_combined)
  # Chord diagram
  pdf(file=filename,
      width = 20, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_chord_gene(object.list[[2]], 
                       slot.name = 'net',
                       sources.use = unique(net.up$source), 
                       targets.use = unique(net.up$target), 
                       net = net.up, 
                       lab.cex = 0.8, small.gap = 0.5, 
                       big.gap = 5,
                       title.name = paste0("Up-regulated ",add_str_ligand_or_receptor," signaling in ", names(object.list)[2]),
                       color.use=colors)
  netVisual_chord_gene(object.list[[1]], 
                       slot.name = 'net',
                       sources.use = unique(net.down$source), 
                       targets.use = unique(net.down$target), 
                       net = net.down, 
                       lab.cex = 0.8, small.gap = 0.5, 
                       big.gap = 5,
                       title.name = paste0("Down-regulated ",add_str_ligand_or_receptor," signaling in ", names(object.list)[2]),
                       color.use=colors)
  #netVisual_chord_gene(object.list[[2]], sources.use = levels(cellchat_combined@idents$joint), targets.use = levels(cellchat_combined@idents$joint), slot.name = 'net', net = net.up, lab.cex = 0.5, small.gap = 2.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),color.use=colors)
  #netVisual_chord_gene(object.list[[1]], sources.use = levels(cellchat_combined@idents$joint), targets.use = levels(cellchat_combined@idents$joint), slot.name = 'net', net = net.down, lab.cex = 0.5, small.gap = 2.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]),color.use=colors)
  dev.off()
  
}


visualize_ranked_signaling_pathways_between_CT_pairs <- function(cellchat_combined,results_path,source_celltypes,target_celltypes){
  #rank signalling changes by information flow or number of interactions
  for (i in seq(1,length(source_celltypes))){
    visualize_ranked_signaling_pathways(cellchat_combined, "weight", source_celltypes[i], target_celltypes[i], results_path)
    visualize_ranked_signaling_pathways(cellchat_combined, "count", source_celltypes[i], target_celltypes[i], results_path)
  }
}



netClustering_multisession <- function (object, slot.name = "netP", type = c("functional", 
                                               "structural"), comparison = NULL, k = NULL, methods = "kmeans", 
          do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, 
          k.eigen = NULL) 
{
  type <- match.arg(type)
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Classification learning of the signaling networks for a single dataset", 
        "\n")
  }
  else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Classification learning of the signaling networks for datasets", 
        as.character(comparison), "\n")
  }
  comparison.name <- paste(comparison, collapse = "-")
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  data.use <- Y
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use, k, nstart = 10)$cluster
    }
    else {
      N <- nrow(data.use)
      kRange <- seq(2, min(N - 1, 10), by = 1)
      if (do.parallel) {
        future::plan("multisession", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(test = future::nbrOfWorkers() == 
                            1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
      results = my.sapply(X = 1:length(kRange), FUN = function(x) {
        idents <- kmeans(data.use, kRange[x], nstart = 10)$cluster
        clusIndex <- idents
        adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, 
                                                   clusIndex, FUN = "==")), nrow = N, ncol = N)
        return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
      }, simplify = FALSE)
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce("+", adjMat)/length(kRange)
      res <- computeEigengap(as.matrix(CM))
      numCluster <- res$upper_bound
      clusters = kmeans(data.use, numCluster, nstart = 10)$cluster
      if (do.plot) {
        gg <- res$gg.obj
        ggsave(filename = paste0("estimationNumCluster_", 
                                 fig.id, "_", type, "_dataset_", comparison.name, 
                                 ".pdf"), plot = gg, width = 3.5, height = 3, 
               units = "in", dpi = 300)
      }
    }
  }
  else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D) - A
    L <- diag(D^-0.5) %*% L %*% diag(D^-0.5)
    evL <- eigen(L, symmetric = TRUE)
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[, (ncol(evL$vectors) - k.eigen + 1):ncol(evL$vectors)]
    clusters = kmeans(Z, k, nstart = 20)$cluster
  }
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
    methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
  return(object)
}

xxx<- function(cellchat,df.net,results_path){
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net") 
  communication_filtered = df.net[order(df.net$prob, decreasing=TRUE),]
  communication_filtered = communication_filtered[communication_filtered$prob > 0.03,] # prioritize by probability of communication
  
  #build color palette for cell types
  #to do make work for our data set
  humouse = read.table("mouse_human_2.txt", sep="\t")
  colnames(humouse) = c("mouse", "human")
  custom_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(2, "Accent"))
  color_count <- nrow(humouse)
  if(color_count > length(custom_palette)) {
    stop("Not enough unique colors in the custom palette for each data row!")
  } else {
    humouse$color <- custom_palette[1:color_count]
  }
  
  cn = colnames(humouse)
  colnames(humouse) = cn
  humouse
}

get_pathways <- function(cellchat_SCZ,cellchat_CTRL){
  df_signalingpws_SCZ <- subsetCommunication(cellchat_SCZ,slot.name = "netP" )
  signalingpws_SCZ  <- unique(df_signalingpws_SCZ$pathway_name)
  df_signalingpws_CTRL <- subsetCommunication(cellchat_CTRL,slot.name = "netP" )
  signalingpws_CTRL  <- unique(df_signalingpws_CTRL$pathway_name)
  pathways.show <- unique(c(signalingpws_SCZ,signalingpws_CTRL))
  return(pathways.show)
}

get_ct_colors <- function(cellchat_combined, with_white=FALSE){
  CTs <- levels(cellchat_combined@idents$joint)
  #[1] "Astrocytes"         "Endothelial/mural"  "Exc L2-3 IT I"      "Exc L2-3 IT II"    
  #[5] "Exc L3-4 IT"        "Exc L3-6 IT"        "Exc L5-6 CT and NP" "Exc L5-6 IT I"     
  #[9] "Exc L5-6 IT II"     "Inh LAMP5"          "Inh PVALB"          "Inh SST"           
  #[13] "Inh VIP"            "Microglia"          "Oligodendrocytes"   "OPCs" 
  colors <- c("#665C47","#66341E","#2EBF5E","#21EC1D","#5959AD","#58D2CF","#0D5A8B","#9AB7E2","#394FD3","#DA808C","#D93137","#FF9900","#B864CC","#94AF97","#53776C","#697255")
  names(colors) <- CTs
  if (with_white==TRUE){
    colors <- c(colors,"#FFFFFF")
    names(colors) <- c(CTs,"white")
  }
  return(colors)
}

get_top5_interactions <- function(df, opt_ligand, opt_up){
  if (opt_ligand==TRUE){
    if (opt_up==TRUE){
      top_interactions <- df %>%
        arrange(desc(ligand.logFC)) %>%  # ascending order
        distinct(interaction_name, .keep_all = TRUE) %>%  # keep first occurrence of each interaction_name
        slice_head(n = 5) %>%  # take top 5
        select(interaction_name) 
    } else{
      top_interactions <- df %>%
        arrange(ligand.logFC) %>%  # ascending order
        distinct(interaction_name, .keep_all = TRUE) %>%  # keep first occurrence of each interaction_name
        slice_head(n = 5) %>%  # take top 5
        select(interaction_name)  # keep only the column you want
    }
  } else{
    if (opt_up==TRUE){
      top_interactions <- df %>%
        arrange(desc(receptor.logFC)) %>%  # ascending order
        distinct(interaction_name, .keep_all = TRUE) %>%  # keep first occurrence of each interaction_name
        slice_head(n = 5) %>%  # take top 5
        select(interaction_name) 
    } else{
      top_interactions <- df %>%
        arrange(receptor.logFC) %>%  # ascending order
        distinct(interaction_name, .keep_all = TRUE) %>%  # keep first occurrence of each interaction_name
        slice_head(n = 5) %>%  # take top 5
        select(interaction_name)  # keep only the column you want
    }
  }
  print(top_interactions)
  return(top_interactions)
}

save_sign_ligand_and_target_interactions_T6 <- function(net.ligand,net.receptor,results_path){

  # Create a new workbook
  wb <- createWorkbook()
  
  # Add sheets and write data
  addWorksheet(wb, "sign_ligands")
  writeData(wb, sheet = "sign_ligands", net.ligand, rowNames = FALSE)
  
  # You can reuse the same dataframes with different filters or summaries
  addWorksheet(wb, "sign_receptors")
  writeData(wb, sheet = "sign_receptors", net.receptor, rowNames = FALSE)
  
  # Save the workbook
  saveWorkbook(wb, file = paste0(results_path,"T6_sign_ligands_and_receptors.xlsx"), overwrite = TRUE) 
}