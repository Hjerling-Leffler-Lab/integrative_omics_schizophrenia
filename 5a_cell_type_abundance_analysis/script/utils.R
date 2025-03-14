getPValueDf <- function(cao, cell.group.order) {
  
  freqs <- cao$test.results$coda$cnts %>% {. / rowSums(.)}
  pval.df <- cao$sample.groups %>% {split(names(.), .)} %>% 
    {matrixTests::col_wilcoxon_twosample(freqs[.[[1]],], freqs[.[[2]],])} %$% 
    setNames(pvalue, rownames(.)) %>% 
    p.adjust("BH")
  #save pval.df
  write.csv(pval.df, file=paste0(results_path, "Pval_df_wilcoxon_freq.csv"))
  pval.df.stars <- pval.df %>% cacoa:::pvalueToCode(ns.symbol="") %>% 
    {tibble(ind=names(.), freq=., coda=cao$test.results$coda$padj[names(.)])} %>% 
    mutate(ind=factor(ind, levels=cell.group.order), coda=cacoa:::pvalueToCode(coda, ns.symbol="")) %>% 
    rename(Freqs=freq, CoDA=coda)
  return(pval.df.stars)
}

addPvalueToCoda <- function(gg, cao, x.vals, show.legend=FALSE, size=4, legend.title="Significance", results_path) {
  pval.df.stars <- getPValueDf(cao, cell.group.order=levels(gg$data$ind))
  #save pval.df
  pval.df_2 <- cao$test.results$coda$padj
  write.csv(pval.df_2, file=paste0(results_path, "Pval_adj_df_2.csv"))
  #save pval.df
  #browser()
  pval.df_2 <- cao$test.results$coda$pval
  write.csv(pval.df_2, file=paste0(results_path, "Pval_df_2.csv"))
  
  gg <- gg + 
    geom_text(aes(x=x.vals[1], label=CoDA, color="CoDA"), data=pval.df.stars, vjust=0.75, size=size) #+
    #geom_text(aes(x=x.vals[2], label=Freqs, color="Wilcox"), data=pval.df.stars, vjust=0.75, size=size) +
    #scale_color_manual(values=c("black", "darkred"))
  
  if (show.legend) {
    gg <- gg + 
      cacoa:::theme_legend_position(c(1, 0.04)) +
      guides(fill=guide_none(), color=guide_legend(title=legend.title))
  }
  return(gg)
}

get_color_palettes <- function(data_path){
  palettes <- c()
  palettes$scz_ctrl_pal <- c(SCZ="#FF8C00" , CTRL="#4682B4")
  
  #load colors for cell types from sheet
  colors_CTs = read.xlsx(paste0(data_path,'Cell_type_colors.xlsx'), sheet='Sheet1')
  palettes$ct_palette = eval(parse(text=paste0("colors_CTs$`Cell.type.(",as.character(n_cluster),").color`")))
  names(palettes$ct_palette) = eval(parse(text=paste0("colors_CTs$`Cell.type.(",as.character(n_cluster),")`")))
  
  return(palettes)
}

#load data from loom
load_data_for_cacoa <- function(data_path, filename){
  
  file_with_path <- paste0(data_path,filename)
  print(file_with_path)
  #Connect Loom File
  #setwd(data_path)
  #data_HC <- connect(file,mode = "r", skip.validate = TRUE)
  data <- connect(file_with_path,mode = "r", skip.validate = TRUE)
  
  return(data)
}

build_cacoa_object <- function(data,n_cluster,n_cores,opt_cluster_based, CT_i, palettes, opt_run_per){
  
  cell_groups <- data[[paste0(paste0('col_attrs/cluster_name_',n_cluster),'CTs')]][]
  #remove CT none (removed)
  if (opt_run_per!='celltype'){
    idx <- which(!(cell_groups %in% c('none (removed)')))
  } else if (opt_run_per == 'celltype'){
    #get CT list
    CT_list = unique(cell_groups)
    if (opt_cluster_based==FALSE){
      #only store values for current CT in object
      idx <- which(cell_groups == CT_list[CT_i])
    } else{
      #idx <- seq(1,length(cell_groups))
      idx <- which(!(cell_groups %in% c('none (removed)')))
    }
  }
  cell_groups <- cell_groups[idx]
  cellIDs <- data[['col_attrs/CellID']][idx]
  sample_groups <- data[['col_attrs/Disease']][idx]
  #sample_groups <- data[['col_attrs/Sex']][]
  
  sample_per_cell <- data[['col_attrs/Donor']][idx]
  names(sample_groups) <- sample_per_cell
  names(cell_groups) <- cellIDs
  names(sample_per_cell) <- cellIDs
  
  genes <- data[['row_attrs/Gene']][]
  #data matrix selection:
  mat <- t(data[["matrix"]][idx,])
  
  rownames(mat) <- genes
  colnames(mat) <- cellIDs
  # make Seurat oject
  seurat_obj <- CreateSeuratObject(counts = mat, project = "SeuratProject", assay="RNA", min.cells = 0, min.features = 0,meta.data = NULL)
  #DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj<-FindVariableFeatures(seurat_obj,nfeatures=2000)
  seurat_obj <- ScaleData(seurat_obj, verbose=FALSE)
  seurat_obj <- RunPCA(seurat_obj,npcs=30,verbose=FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims=1:25)
  #how much var is explained?
  #(seurat_obj@reductions$pca@stdev)^2
  if (opt_cluster_based==FALSE){
    seurat_obj<- FindNeighbors(seurat_obj, reduction="pca", dims = 1:25, k.param = 60, prune.SNN = 1/15) #reduction= "umap", dims = 1:2, verbose = T)
    seurat_obj<- FindClusters(seurat_obj, algorithm= 1, resolution = 0.5, verbose = T, graph.name = "RNA_snn")
  }
  
  #will be automatically stored in cao object
  
  #build the cacoa object
  cao <- cacoa::Cacoa$new(
    seurat_obj, sample.groups=sample_groups, cell.groups=cell_groups,sample.per.cell = sample_per_cell,
    target.level='SCZ', ref.level='CTRL', n.cores=n_cores, verbose=FALSE, graph.name="RNA_snn"
  )
  rm(seurat_obj)
  
  #store colors:
  cao$sample.groups.palette = palettes$scz_ctrl_pal
  cao$cell.groups.palette = palettes$ct_palette
  
  #cao$embedding <- Embeddings(seurat_obj,reduction="pca")[,1:20]
  return(cao)
}

#resort CT labels
resort_CT_labels <- function(CT_levels){
  CT_levels_ordered <- CT_levels[CT_levels!='none (removed)'] 
  idx_micro <- which(CT_levels_ordered=="Microglial cells")
  idx_endo <- which(CT_levels_ordered=="Endothelial and mural cells" )
  CT_levels_ordered <- CT_levels_ordered[c(seq(idx_endo+1,idx_micro-1),seq(1,idx_endo),seq(idx_micro,length(CT_levels_ordered)))]
  return(CT_levels_ordered)
}

cluster_based_analysis <- function(cao, results_path){
  #resort CT labels 
  if (opt_inhibitory_only==FALSE){
    CT_levels_ordered <- resort_CT_labels(levels(cao$cell.groups))
  } else{
    CT_levels_ordered <- levels(cao$cell.groups)[levels(cao$cell.groups)!='none (removed)'] 
  }
  
  graphic_filename_1 <- paste0(results_path,"clusterBased_group_size.pdf")
  pdf(file = graphic_filename_1,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(cao$plotCellGroupSizes(show.significance=TRUE, legend.position='top')+coord_flip()+scale_x_discrete(limits = rev(CT_levels_ordered)))
  dev.off()
  
  #plot compositional changes:
  cao$estimateCellLoadings()
  p_val_per_CT <- cao$test.results$coda$pval
  pval_smallest <- as.numeric(p_val_per_CT[p_val_per_CT==min(p_val_per_CT)])
  update_geom_defaults("text", list(size = 3))
  
  gg_loadings <- cao$plotCellLoadings(show.pvals=FALSE, alpha=0.1, annotation.x=1.0) + scale_x_continuous(limits=c(-1.0,1.2), expand=c(0, 0.0, 0.0, 0.1)) #+xlim(c(-1,1)
  gg_loadings %<>% addPvalueToCoda(cao, c(0.53, 0.7), show.legend=TRUE, size=4, legend.title="Significance", results_path=results_path)
  graphic_filename_2 <- paste0(results_path,"clusterBased_loadings.pdf")
  pdf(file = graphic_filename_2,   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(gg_loadings) #+ geom_text(x=1, y=length(CT_levels_ordered), label=paste0("p-value = ",pval_smallest))) #The red line here shows statistical significance
  dev.off()
  # if there is no red line --> no significant changes
  # 
  #plot contrast tree
  tree_theme <- theme(
  legend.key.height=unit(10, "pt"), legend.key.width=unit(14, "pt"),
  legend.position="bottom", plot.margin=margin(),
  axis.text.y=element_text(hjust=1, vjust=0.5, margin=margin()), axis.text.x=element_blank(),
  axis.ticks=element_blank()
  )
  #browser()
  counts_per_sample_and_CT <- cao$test.results$coda$cnts
  group_per_sample <- cao$sample.groups[!duplicated(names(cao$sample.groups))]
  sample_names <- names(group_per_sample)
  group_per_sample <- unfactor(group_per_sample)
  names(group_per_sample) <- sample_names
  gg_ms_tree <- plotContrastTree_v2(d.counts=counts_per_sample_and_CT, d.groups=group_per_sample, target.level='SCZ', ref.level='CTRL')

  gg_ms_tree <- cao$plotContrastTree() + coord_flip() + tree_theme + guides(color=guide_legend(direction="vertical", title="Condition", order=1), fill=guide_colorbar(title.position="top", order=2, title.hjust=0.5))
  graphic_filename_3 <- paste0(results_path,"clusterBased_loadings_contrast_tree.pdf")
  pdf(file = graphic_filename_3,   # The directory you want to save the file in
      width = 8, # The width of the plot in inches
      height = 8) # The height of the plot in inches
  print(gg_ms_tree)
  dev.off()
  
  #plot expression changes:
  cao$estimateExpressionShiftMagnitudes(n.cores=n_cores)
  
  graphic_filename_4 <- paste0(results_path,"clusterBased_expr_shift_magnitudes.pdf")
  pdf(file = graphic_filename_4,   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(cao$plotExpressionShiftMagnitudes()) #Here, y-axis shows magnitude of changes, while asterisks on top of bars show their significance.
  dev.off()
  #Here, y-axis shows magnitude of changes, while asterisks on top of bars show their significance.
  #analyse influence of metadata (only meaningful if something is significant)
}

plotContrastTree_v2 <- function (d.counts, d.groups, ref.level, target.level, plot.theme = theme_get(), 
          label.angle = 90, p.threshold = 0.05, adjust.pvalues = TRUE, 
          h.methods = "both", font.size = 3, label.hjust = 1, show.text = "pvalue.code", 
          tree.order = NULL, loadings.mean = NULL, palette = NULL, 
          verbose = FALSE) 
{
  #checkPackageInstalled(c("ggdendro"), cran = TRUE)
  log.f <- getLogFreq(d.counts)
  if (h.methods == "up") {
    if (verbose) 
      message("up")
    t.cur <- constructTreeUp(d.counts, d.groups)
  }
  else if (h.methods == "down") {
    if (verbose) 
      message("down")
    t.cur <- constructTree(d.counts, d.groups)
  }
  else {
    if (verbose) 
      message("up and down")
    t.cur <- constructTreeUpDown(d.counts, d.groups)
  }
  if (!is.null(loadings.mean) && is.null(tree.order)) {
    tree.order <- names(sort(loadings.mean))
  }
  if (!is.null(tree.order)) {
    t <- t.cur$tree
    tree.order <- intersect(tree.order, t$tip.label)
    d <- distTreeOrder(t, tree.order)
    idx <- min(t$edge[, 1]):max(t$edge[, 1])
    for (i.node in idx) {
      t.alt <- ape::rotate(t, i.node)
      d.alt <- distTreeOrder(t.alt, tree.order)
      if (d.alt <= d) {
        t <- t.alt
        d <- d.alt
      }
    }
    t.cur$tree <- t
    t.cur$dendro <- compute.brlen(t.cur$tree, method = "Grafen") %>% 
      as.hclust() %>% as.dendrogram()
  }
  tree <- t.cur$tree
  sbp <- sbpInNodes(tree)
  dend.data <- ggdendro::dendro_data(t.cur$dendro, type = "rectangle")
  types.order <- dend.data$labels$label
  for (k in 1:ncol(sbp)) {
    p <- sbp[, k]
    type.plus <- rownames(sbp)[p > 0]
    type.minus <- rownames(sbp)[p < 0]
    if (which(types.order == type.plus[1]) < which(types.order == 
                                                   type.minus[1])) {
      sbp[, k] <- -sbp[, k]
    }
  }
  node.pos <- dend.data$segments %$% .[(y == yend) & (yend != 
                                                        0), ]
  node.pos$id <- tree$edge[, 1]
  node.pos$to <- tree$edge[, 2]
  innode.pos <- unique(node.pos[, c("x", "y", "id")])
  rownames(innode.pos) <- innode.pos$id
  innode.pos$range <- -1
  for (i in 1:nrow(innode.pos)) {
    tmp <- node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] <- max(tmp) - min(tmp)
  }
  innode.pos <- innode.pos[order(innode.pos$id), ]
  balances <- getNodeBalances(log.f, sbp)
  colnames(balances) <- rownames(innode.pos)
  p.val <- c()
  for (i in 1:ncol(balances)) {
    aov.data <- data.frame(balance = balances[, i], group = d.groups)
    mod <- lm(group ~ balance, data = aov.data)
    res <- summary(mod)
    p.val <- c(p.val, res$coefficients[2, 4])
  }
  p.val[is.na(p.val)] <- 1
  p.adj <- if (adjust.pvalues) 
    p.adjust(p.val, method = "fdr")
  else p.val
  px.init <- createDendrogram(dend.data, plot.theme = plot.theme, 
                              font.size = font.size, angle = label.angle, hjust = label.hjust)
  if (sum(p.adj < p.threshold) == 0) 
    return(px.init)
  df.pval <- data.frame()
  df.bals <- data.frame()
  df.bal.quartile <- data.frame()
  df.bal.median <- data.frame()
  df.text <- data.frame()
  group.levels <- c(ref.level, target.level)
  for (id.node in 1:ncol(balances)) {
    if (p.adj[id.node] > p.threshold) 
      next
    df.pval <- rbind(df.pval, innode.pos[id.node, c("x", 
                                                    "y")])
    x.tmp <- balances[, id.node]
    x.tmp <- x.tmp - mean(x.tmp)
    df.text %<>% rbind(data.frame(x = innode.pos$x[id.node] + 
                                    innode.pos$range[id.node]/2, y = innode.pos$y[id.node], 
                                  range = sprintf("%2.1f", max(abs(x.tmp))), pvalue = signif(p.adj[id.node], 
                                                                                             2), pvalue.code = pvalueToCode(p.adj[id.node])))
    x.tmp <- x.tmp/max(abs(x.tmp))/2 * innode.pos$range[id.node] * 
      0.9
    x.tmp <- x.tmp + innode.pos$x[id.node]
    y.tmp <- d.groups * 0.03 + innode.pos$y[id.node] - 0.05
    q.case <- quantile(x.tmp[d.groups], c(0.25, 0.75))
    q.control <- quantile(x.tmp[!d.groups], c(0.25, 0.75))
    df.bal.quartile %<>% rbind(data.frame(x = c(q.case, 
                                                q.control), y = c(rep(y.tmp[d.groups][1], 2), rep(y.tmp[!d.groups][1], 
                                                                                                  2)), group = group.levels[1 + c(1, 1, 0, 0)], node = id.node))
    df.bal.median %<>% rbind(data.frame(x = c(median(x.tmp[d.groups]), 
                                              median(x.tmp[!d.groups])), y = c(y.tmp[d.groups][1], 
                                                                               y.tmp[!d.groups][1]), group = group.levels[1 + c(1, 
                                                                                                                                0)], node = id.node))
    df.bals %<>% rbind(data.frame(x = x.tmp, y = y.tmp, 
                                  group = group.levels[1 + d.groups], node = id.node))
  }
  px <- px.init + geom_point(data = df.bals, aes(x = x, y = y, 
                                                 col = as.factor(group), group = as.factor(node)), alpha = 0.1, 
                             size = 1) + geom_point(data = df.bal.median, aes(x = x, 
                                                                              y = y, col = as.factor(group)), size = 2.5, shape = 18) + 
    geom_point(data = df.pval, aes(x = x, y = y)) + geom_line(data = df.bal.quartile, 
                                                              aes(x = x, y = y, col = as.factor(group), group = interaction(group, 
                                                                                                                            node)), size = 0.75) + labs(col = " ")
  if (is.logical(show.text) && !show.text) 
    show.text <- NULL
  if (!is.null(show.text)) {
    if (!(show.text %in% colnames(df.text))) 
      stop("Unexpected value for show.text: ", show.text)
    px <- px + geom_text(data = df.text, mapping = aes_string(x = "x", 
                                                              y = "y", label = show.text), vjust = 0, hjust = 0, 
                         size = font.size)
  }
  if (!is.null(loadings.mean)) {
    node.leaves <- node.pos[node.pos$to < min(tree$edge[, 
                                                        1]), ]
    node.leaves$label <- tree$tip.label[node.leaves$to]
    node.leaves$loadings <- loadings.mean[node.leaves$label]
    px <- px + geom_tile(data = node.leaves, mapping = aes(x = xend, 
                                                           y = -0.1, fill = loadings, width = 0.5, height = 0.1)) + 
      guides(fill = guide_colorbar(title = "loadings", 
                                   title.position = "top", direction = "horizontal", 
                                   title.hjust = 0.5))
  }
  if (!is.null(palette)) {
    px <- px + scale_color_manual(values = palette) + scale_fill_gradient2(low = palette[ref.level], 
                                                                           high = palette[target.level], mid = "grey80", midpoint = 0)
  }
  return(px)
}

getLogFreq <- function (d.counts) 
{
  checkData(d.counts)
  data.norm <- d.counts
  data.norm[data.norm == 0] <- 1
  log.freq <- t(apply(data.norm, 1, function(y) log(y/sum(y))))
  return(log.freq)
}

checkData <- function (d.counts) 
{
  if (is.null(d.counts)) 
    stop("Cell count matrix is not provided")
  if (nrow(d.counts) == 0) 
    stop("Sample size is zero")
  if (ncol(d.counts) == 0) 
    stop("Cell count information is missed")
}

constructTreeUpDown <- function (cnts, groups) 
{
  cnts[cnts <= 0] <- 0.5
  ref.set <- referenceSet(cnts, groups)
  cell.lists <- ref.set$cell.list
  cell.types <- colnames(cnts)
  freqs <- (cnts)/rowSums(cnts)
  freqs.lists <- c()
  for (i in 1:length(cell.lists)) {
    freqs.lists <- cbind(freqs.lists, apply(freqs[, cell.lists[[i]], 
                                                  drop = FALSE], 1, psych::geometric.mean))
  }
  colnames(freqs.lists) <- paste("tmp", 1:length(cell.lists), 
                                 sep = "")
  t.list <- constructBestPartitionTree(freqs.lists, groups)
  sbp.list <- t.list$sbp
  sbp.all <- matrix(ncol = 0, nrow = length(cell.types), dimnames = list(cell.types, 
                                                                         c()))
  for (k in 1:ncol(sbp.list)) {
    p <- sbp.list[, k]
    i.list.plus <- which(p > 0)
    i.list.minus <- which(p < 0)
    sbp.tmp <- rep(0, length(cell.types))
    for (i.plus in i.list.plus) {
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = 1
    }
    for (i.plus in i.list.minus) {
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = -1
    }
    sbp.all <- cbind(sbp.all, sbp.tmp)
  }
  for (i in 1:length(cell.lists)) {
    freqs.tmp <- freqs[, cell.lists[[i]], drop = FALSE]
    if (ncol(freqs.tmp) == 1) 
      next
    tmp <- constructBestPartitionTree(freqs.tmp, groups)
    sbp.tmp <- tmp$sbp
    sbp.add <- matrix(0, nrow = nrow(sbp.all), ncol = ncol(sbp.tmp))
    rownames(sbp.add) <- rownames(sbp.all)
    sbp.add[rownames(sbp.tmp), ] <- sbp.tmp
    sbp.all <- cbind(sbp.all, sbp.add)
  }
  tree <- sbp2tree(sbp.all)
  h.tmp <- compute.brlen(tree, method = "Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  return(list(tree = tree, sbp = sbp.all, dendro = d.cur))
}

referenceSet <- function (freqs, groups, p.thresh = 0.05) 
{
  checkPackageInstalled("psych", cran = TRUE)
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)
  for (it in 1:(length(cell.list) - 2)) {
    mx <- matrix(0, nrow = length(cell.list), ncol = length(cell.list))
    for (i in 1:length(cell.list)) {
      for (j in 1:length(cell.list)) {
        if (j <= i) 
          next
        ratio <- log(apply(freqs[, cell.list[[i]], drop = FALSE], 
                           1, psych::geometric.mean)) - log(apply(freqs[, 
                                                                        cell.list[[j]], drop = FALSE], 1, psych::geometric.mean))
        
        #to do: add ratio and groups to a dataframe, using their names as matching variable
        df = data.frame(groups)
        df$groups <- as.factor(df$groups)
        df$ratio <- ratio[match(row.names(df),names(ratio))]
        #make sure not nans or infs exists
        #Replace NaN & Inf with NA
        df[is.na(df) | df=="Inf"] = NA
        mod <- lm_v2(ratio ~ groups, data = df)
        res <- summary(mod)
        pval <- res$coefficients[2, 4]
        mx[i, j] <- pval
      }
    }
    
    if (max(max(mx)) < p.thresh) 
      break
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], 
                                            cell.list[[j]])
    cell.list <- cell.list[-c(i, j)]
    if (it == 1) {
      mx.first <- mx
      rownames(mx.first) <- cell.types
      colnames(mx.first) <- cell.types
    }
  }
  len.cell.list <- sapply(cell.list, length)
  id.ref <- which(len.cell.list == max(len.cell.list))
  for (i in 1:nrow(mx.first)) {
    for (j in 1:ncol(mx.first)) {
      if (j <= i) 
        next
      mx.first[j, i] <- mx.first[i, j]
    }
  }
  return(list(ref.set = cell.list[[id.ref[1]]], cell.list = cell.list, 
              mx.first = mx.first, mx.final = mx))
}

lm_v2 <- function(formula, data, subset, weights, na.action, method = "qr", 
                            model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
                            contrasts = NULL, offset, ...) 
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame") 
    return(mf)
  else if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response_v2(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  offset <- model.offset(mf)
  mlm <- is.matrix(y)
  ny <- if (mlm) 
    nrow(y)
  else length(y)
  if (!is.null(offset)) {
    if (!mlm) 
      offset <- as.vector(offset)
    if (NROW(offset) != ny) 
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    NROW(offset), ny), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (mlm) matrix(NA_real_, 0, 
                                             ncol(y)) else numeric(), residuals = y, fitted.values = 0 * 
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
                                                                                0) else ny)
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    #print(y)
    #print(x)
    z <- if (is.null(w)) 
      lm.fit_v2(x, y, offset = offset, singular.ok = singular.ok, 
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
                 ...)
  }
  class(z) <- c(if (mlm) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model) 
    z$model <- mf
  if (ret.x) 
    z$x <- x
  if (ret.y) 
    z$y <- y
  if (!qr) 
    z$qr <- NULL
  z
}

model.response_v2 <- function (data, type = "any") 
{
  if (attr(attr(data, "terms"), "response")) {
    if (is.list(data) | is.data.frame(data)) {
      v <- data[[1L]]
      if (type == "numeric" && is.factor(v)) {
        warning("using type = \"numeric\" with a factor response will be ignored")
      }
      else if (type == "numeric" | type == "double") 
        storage.mode(v) <- "double"
      else if (type != "any") 
        stop("invalid response type")
      if (is.matrix(v) && ncol(v) == 1L) 
        dim(v) <- NULL
      rows <- attr(data, "row.names")
      if (nrows <- length(rows)) {
        if (length(v) == nrows) 
          names(v) <- rows
        else if (length(dd <- dim(v)) == 2L) 
          if (dd[1L] == nrows && !length((dn <- dimnames(v))[[1L]])) 
            dimnames(v) <- list(rows, dn[[2L]])
      }
      return(v)
    }
    else stop("invalid 'data' argument")
  }
  else return(NULL)
}

checkPackageInstalled <- function (pkgs, details = "to run this function", install.help = NULL,
          bioc = FALSE, cran = FALSE)
{
  pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(pkgs) == 0) {
    return(NULL)
  }
  if (length(pkgs) > 1) {
    pkgs <- paste0("c('", paste0(pkgs, collapse = "', '"),
                   "')")
    error.text <- paste("Packages", pkgs, "must be installed",
                        details)
  }
  else {
    pkgs <- paste0("'", pkgs, "'")
    error.text <- paste(pkgs, "package must be installed",
                        details)
  }
  if (!is.null(install.help)) {
    error.text <- paste0(error.text, ". Please, run `",
                         install.help, "` to do it.")
  }
  else if (bioc) {
    error.text <- paste0(error.text, ". Please, run `BiocManager::install(",
                         pkgs, ")", "` to do it.")
  }
  else if (cran) {
    error.text <- paste0(error.text, ". Please, run `install.packages(",
                         pkgs, ")", "` to do it.")
  }
  stop(error.text)
}

#' Construct the canonical tree
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructBestPartitionTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
  # ---------

  n.cells <- ncol(cnts)

  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(cnts), c()))

  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list for separation contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    # ------------------------------------------------
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]

    # -------
    # Find the best partition
    sbp.best <- bestPartition(d.tmp, groups)

    # -------
    # Get cell types from opposite sides of the principal balance
    cells.tmp <- colnames(d.tmp)

    cells.plus <- cells.tmp[sbp.best > 0]
    cells.minus <- cells.tmp[sbp.best < 0]
    #
    # print(cells.plus)
    # print(cells.minus)

    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1

    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }

    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }

  }

  # res <- list(sbp = sbp.cda, partiotions = unsolved.cells)
  tree <- sbp2tree(sbp.cda)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)

  return(list(tree = tree, sbp = sbp.cda, dendro = d.cur))
}

#' Check groups
#'
#' @param d.groups Groups variable for samples
#' @keywords internal
checkGroups <- function(d.groups){
  if(is.null(d.groups))
    stop('Groups are not provided')
  if(length(d.groups) == 0)
    stop('Group information is missed')
}

#' Check cell count data and groups
#'
#' @param d.counts Cell count table
#' @param d.groups Groups variable for samples
#' @keywords internal
checkDataGroups <- function(d.counts, d.groups){
  checkData(d.counts)
  checkGroups(d.groups)
  if(nrow(d.counts) != length(d.groups))
    stop('Sample size in cell count matrix and in group variable is not the same')
}


#' Get tree from balances
#'
#' @param sbpart A contrast matrix for balances, sequential binary partition
#' @return phylo tree
#' @details We use ellipsis because a previous version of this function
#'          contains additional and insignificant parameters.
#' @keywords internal
sbp2tree <- function(sbpart){
  checkSbpWhole(sbpart)
  # ---------

  n.cells <- nrow(sbpart)
  edges <- c(n.cells+1, n.cells+1)
  id <- 2*n.cells - 1

  collapse <- sbpart
  rownames(collapse) <- 1:n.cells

  for(i in rev(1:(n.cells - 1)) ){
    ids <- which(collapse[,i] != 0)

    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[1]])))
    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[2]])))
    rownames(collapse)[ids[2]] <- id
    id <- id - 1

    collapse[ids[1],] <- 0
  }

  # Create face tree, because we need a tidytree object
  x <- tidytree::as_tibble(rtree(n = n.cells, br = NULL))
  x[,1:2] <- edges

  t.tmp <- as.phylo(x)
  t.tmp$tip.label <- rownames(sbpart)
  return(t.tmp)
}

#' Check sequential binary partitions(sbp) for the whole set of balances
#'
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbpWhole <- function(sbp){
  checkSbp(sbp)
  if(nrow(sbp) < 2)
    stop('Balances require at least 2 cell types')
  if((nrow(sbp) - ncol(sbp)) != 1)
    stop('Wrong number of balances provided')
}

#' Check the agreement between sbp and cell count table
#'
#' @param d.counts Cell count table
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkDataSbp <- function(d.counts, sbp){
  checkData(d.counts)
  checkSbp(sbp)
  if(ncol(sbp) != ncol(d.counts))
    stop('Number of cell types in data and sbp should be equal')
}

#' Check sequential binary partitions(sbp)
#'
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbp <- function(sbp){
  if(is.null(sbp))
    stop('Sequential binary partitions are not provided')
  if(nrow(sbp) == 0)
    stop('No partitions are provided')
}

#' @keywords internal
bestPartition <- function(freqs.tmp, groups){


  n.cells <- ncol(freqs.tmp)
  sbp <- 1
  for(i in 2:n.cells){
    sbp <- rbind(cbind(sbp, 1),
                 cbind(sbp, -1))
  }
  sbp <- sbp[-c(1, nrow(sbp)),]
  log.f <- log(freqs.tmp)

  for(i in 1:nrow(sbp)){
    p <- sbp[i,]
    p[p < 0] <- p[p < 0] / sum(p < 0)
    p[p > 0] <- p[p > 0] / sum(p > 0)
    sbp[i,] <- p
  }
  colnames(sbp) <- colnames(freqs.tmp)


  bals <- sbp %*% t(log.f)
  p <- c()
  for(i in 1:nrow(bals)){
    b <- bals[i,]
    mod <- lm(groups ~ b)
    res <- summary(mod)
    pval <- res$coefficients[2,4]
    p <- c(p, pval)
  }

  sbp.best <- (sbp[p == min(p),] > 0) * 2 - 1
  names(sbp.best) <- colnames(freqs.tmp)
  # print(sbp.best)


  return(sbp.best)
}


lm.fit_v2 <- function (x, y, offset = NULL, method = "qr", tol = 1e-07, singular.ok = TRUE, 
          ...) 
{
  #browser()
  if (is.null(n <- nrow(x))) 
    stop("'x' must be a matrix")
  if (n == 0L) 
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L) {
    return(list(coefficients = numeric(), residuals = y, 
                fitted.values = 0 * y, rank = 0, df.residual = length(y)))
  }
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1) 
    y <- drop(y)
  if (!is.null(offset)) 
    y <- y - offset
  if (NROW(y) != n) 
    stop("incompatible dimensions")
  if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
                     method), domain = NA)
  chkDots(...)
  z <- .Call(stats:::C_Cdqrls, x, y, tol, FALSE)
  if (!singular.ok && z$rank < p) 
    stop("singular fit encountered")
  coef <- z$coefficients
  pivot <- z$pivot
  r1 <- seq_len(z$rank)
  dn <- colnames(x)
  if (is.null(dn)) 
    dn <- paste0("x", 1L:p)
  nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  r2 <- if (z$rank < p) 
    (z$rank + 1L):p
  else integer()
  if (is.matrix(y)) {
    coef[r2, ] <- NA
    if (z$pivoted) 
      coef[pivot, ] <- coef
    dimnames(coef) <- list(dn, colnames(y))
    dimnames(z$effects) <- list(nmeffects, colnames(y))
  }
  else {
    coef[r2] <- NA
    if (z$pivoted) 
      coef[pivot] <- coef
    names(coef) <- dn
    names(z$effects) <- nmeffects
  }
  z$coefficients <- coef
  r1 <- y - z$residuals
  if (!is.null(offset)) 
    r1 <- r1 + offset
  if (z$pivoted) 
    colnames(z$qr) <- colnames(x)[z$pivot]
  qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]
  c(z[c("coefficients", "residuals", "effects", "rank")], 
    list(fitted.values = r1, assign = attr(x, "assign"), 
         qr = structure(qr, class = "qr"), df.residual = n - 
           z$rank))
}

plot_graph_based_results <- function(cao,CT_name,results_path){
  
  graphic_filename_0 <- paste0(results_path,"clusterFree_density_",CT_name,"_graph.pdf")
  pdf(file = graphic_filename_0,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(cao$plotCellDensity(add.points=FALSE, show.cell.groups=opt_show_cell_groups, legend.title="Density"))
  dev.off()
  
  graphic_filename_1 <- paste0(results_path,"clusterFree_density_diff_",CT_name,"_graph.pdf")
  pdf(file = graphic_filename_1,   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 5) # The height of the plot in inches
  print(plot_grid(
    gg_embedding, 
    gg_z,
    gg_z_adj, 
    ncol=3
  ))
  dev.off()
}

plot_density_based_results <- function(cao,CT_name,results_path,opt_show_cell_groups,gg_embedding,gg_z_kde,gg_z_adj_kde){
  graphic_filename_2 <- paste0(results_path,"clusterFree_density_",CT_name,"_kde.pdf")
  pdf(file = graphic_filename_2,   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  print(cao$plotCellDensity(name='cell.density.kde',add.points=FALSE, show.cell.groups=opt_show_cell_groups, legend.title="Density"))
  dev.off()
  
  
  graphic_filename_3 <- paste0(results_path,"clusterFree_density_diff_",CT_name,"_kde.pdf")
  pdf(file = graphic_filename_3,   # The directory you want to save the file in
      width = 15, # The width of the plot in inches
      height = 5) # The height of the plot in inches
  print(plot_grid(
    gg_embedding, 
    gg_z_kde,
    gg_z_adj_kde, 
    ncol=3
  ))
  dev.off()
}