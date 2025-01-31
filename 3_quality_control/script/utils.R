#author: Lisa Bast
#date: 24.09.2021
#save dataframe from pagoda as loom file

library(loomR)

save_as_loom_file <- function(data, data_filtered, filepath, filename, sample_names_vec, library_names_vec){
  # identify which indices are not filtered out
  #how many cells in data --> if too many look at chunks of genes to make apply still feasible
  n_cells <- dim(data)[2]
  if (n_cells>17000){
    n_genes <- dim(data)[1]
    n_chunks <- 10
    chunk_size <- floor(n_genes/n_chunks)
    for (i in seq(1,n_chunks)){
      #save each column as chunks of genes as a string
      idx_start <-chunk_size*(i-1)+1
      idx_end <- chunk_size*i
      sA_i <- apply(data[idx_start:idx_end,],2,paste,collapse=' ')
      sB_i <- apply(data_filtered[idx_start:idx_end,],2,paste,collapse=' ')
      #compare strings
      if (i==1){
        kept_in_chunk_i <- (sA_i %in% sB_i)
      }else{
        kept_in_chunk_i <- kept_in_chunk_i & (sA_i %in% sB_i)
      }
    }
    #remaining genes are ignored:
    #if (n_genes+1>chunk_size*i){
    #  sA_i <- apply(data[chunk_size*i:n_genes,],2,paste,collapse=' ')
    #  sB_i <- apply(data_filtered[chunk_size*i:n_genes,],2,paste,collapse=' ')
    #  kept_in_chunk_i <- kept_in_chunk_i & (sA_i %in% sB_i)
    #}
    pagoda_filtered <- !(kept_in_chunk_i)
  }else{
    #save each column as a string
    sA <- apply(data,2,paste,collapse=' ')
    sB <- apply(data_filtered,2,paste,collapse=' ')
    #compare strings
    pagoda_filtered <- !(sA %in% sB)
    #Idx_keep <- which(sA %in% sB)
    #pagoda_filtered <- !(Idx_keep %in% seq(1,dim(data)[2]))
  }
  
  #define a cell ID
  print(length(pagoda_filtered))
  CellIDs <- paste0(seq(0,dim(data)[2]-1), sample_names_vec)
  print(length(CellIDs))
  print(length(sample_names_vec))
  print(length(library_names_vec))
  #define new file name
  new_filename <- paste0(filepath,filename)
  #define row and col attributes
  row_attrs <- list(rownames(data))
  names(row_attrs) <- c("Gene")
  col_attrs <- list(CellIDs, sample_names_vec, library_names_vec, pagoda_filtered)
  names(col_attrs) <- c("CellID","Donor","Library","Pagoda filtered")
  #create new loom file
  lfile <- create(new_filename, data, overwrite = TRUE) #columns (cells) by rows (genes)
  # add attributes
  lfile$add.row.attribute(row_attrs,overwrite=TRUE) 
  lfile$add.col.attribute(col_attrs,overwrite=TRUE)
  #close loom file
  lfile$close_all()
}

plot_correlation_of_metrics <- function(metrics, selected_cols,path_results,fig_str,corr_method,opt_plot_ellipses){
  if (opt_plot_ellipses){
    #calculate correlation of matrix
    #plot correlation matrix as ellipses
    corr <- cor(x = metrics[selected_cols],use = "pairwise.complete.obs" ,method=corr_method)
    png(file=paste0(path_results,"QC_metrics_ellipse_plot_",corr_method,fig_str))
    corrplot(corr, method="ellipse")
    dev.off()
  } else{
    #simple ggpairs plot
    png(file=paste0(path_results,"QC_metrics_ggpairs_plot_",corr_method,fig_str))
    ggpairs(
      metrics,
      columns = selected_cols, 
      legend=1,
      mapping=ggplot2::aes(colour=Group),
      title = paste0(corr_method," correlation of quality metrics"),
      diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"), 
      upper = list(continuous = wrap("cor", size = 4), corMethod = corr_method),
      lower = list(continuous = wrap("smooth", alpha=0.3, size=0.1))
    ) + theme(legend.position = "bottom", axis.text.x = element_text(size=12,angle=-90), axis.text.y = element_text(size=12))
    dev.off()
  }
}

plot_PCA_of_samples <- function(metrics,selected_cols,path_results,opt_plot_expl_var,opt_color_by_contrib,opt_color_by){
  #plot pca of log(x+1) transformed values
  metrics_log <- log(metrics[selected_cols]+array(1,c(nrow(metrics[selected_cols]),ncol(metrics[selected_cols]))))
  rownames(metrics_log) <- metrics$donor_ID_python
  res.pca <- prcomp(metrics_log,scale=TRUE)#variables are scaled to have unit variance before the analysis takes place
  
  #percentage explained variance
  if (opt_plot_expl_var){
    png(file=paste0(path_results,"QC_metrics_PCA_percentage_explained_variance",fig_str))
    fviz_eig(res.pca)
    dev.off()
  }
  
  #individuals contribution pc1 vs pc2
  if (opt_color_by_contrib){
    #contrib: color by contribution
    var <- "contrib"
  }else{
    #cos2 : color by quality of representation
    var <- "cos2"
  }
  
  png(file=paste0(path_results,"QC_metrics_PCA_plot_",var,"_individuals",fig_str))
  fviz_pca_ind(res.pca,
               col.ind = var, # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  ) + scale_color_gradient2(low="lightblue",mid="blue", high="red", midpoint=0.55)
  dev.off()
  
  
  if (opt_color_by=="Group"){
    #color by group
    png(file=paste0(path_results,"QC_metrics_PCA_plot_individuals_colored_by_",opt_color_by,fig_str))
    fviz_pca_ind(res.pca,label = "none",
                 habillage = metrics$Group,
                 addEllipses=TRUE, 
                 ellipse.level=0.95)
    dev.off() 
  } else if(opt_color_by=="Library"){
    png(file=paste0(path_results,"QC_metrics_PCA_plot_individuals_colored_by_",opt_color_by,fig_str))
    fviz_pca_ind(res.pca,label = "none",
                 habillage = metrics$Library,
                 addEllipses=TRUE, 
                 ellipse.level=0.95)
    dev.off() 
  }
  
  #variables contribution pc1 vs pc2
  png(file=paste0(path_results,"QC_metrics_PCA_plot_variable_contribution",fig_str))
  fviz_pca_var(res.pca, 
               col.var="contrib"
  )+ scale_color_gradient2(low="lightblue",mid="blue", high="red", midpoint=20) + theme_minimal()
  dev.off()
  
}


# pearson vs spearman correlation
# histograms for each metric are on diagonal
plot_pearson_vs_spearman_correlation <- function(metrics,selected_cols,path_results,fig_str){
  png(file=paste(path_results,paste("QC_metrics_ggpairs_plot_pearson_vs_spearman",fig_str,sep=''),sep=''))
  ggpairs(
    metrics, 
    columns = selected_cols, 
    legend=1,
    mapping=ggplot2::aes(colour=Group),
    title = "Pearson (lower) vs. Spearman (upper) correlation of quality metrics",
    upper = list(continuous = wrap('cor', size=4, method = "spearman")), 
    lower = list(continuous = wrap("cor", size = 4)),
    diag = list(continuous = "barDiag", discrete = "barDiag", na = "naDiag")
  ) + theme(legend.position = "bottom")
  dev.off()
}
