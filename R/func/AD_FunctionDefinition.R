#####Function Definition: Initialization, plotting and p-Tau analysis#####
#####Initialize#####
ad_ini <- function(){
  suppressMessages(require(Seurat))
  suppressMessages(require(SeuratDisk))
  suppressMessages(require(raster))
  suppressMessages(require(stringr))
  suppressMessages(require(dplyr))
  suppressMessages(require(parallel))
  suppressMessages(require(pbapply))
  suppressMessages(require(ggplot2))
  suppressMessages(require(readr))
  suppressMessages(require(EBImage))
  suppressMessages(require(tiff))
  suppressMessages(require(ggpubr))
  suppressMessages(require(monocle3))
  suppressMessages(library(factoextra))
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(umap))
  print("Initialization Finished")
}


#####General#####  
  
ad_plot_cells <- function (cds, x = 1, y = 2, 
                           reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
                           color_cells_by = "cluster", group_cells_by = c("cluster", "partition"), 
                           genes = NULL, show_trajectory_graph = TRUE, 
          trajectory_graph_color = "grey28", trajectory_graph_segment_size = 0.75, 
          norm_method = c("log", "size_only"), label_cell_groups = TRUE, 
          label_groups_by_cluster = TRUE, group_label_size = 2, labels_per_group = 1, 
          label_branch_points = TRUE, label_roots = TRUE, label_leaves = TRUE, 
          graph_label_size = 2, cell_size = 0.35, cell_stroke = I(cell_size/2), 
          alpha = 1, min_expr = 0.1, rasterize = FALSE, scale_to_range = FALSE) 
{
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]), 
                          msg = paste("No dimensionality reduction for", reduction_method, 
                                      "calculated.", "Please run reduce_dimensions with", 
                                      "reduction_method =", reduction_method, "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >= max(x, y), 
                          msg = paste("x and/or y is too large. x and y must", 
                                      "be dimensions in reduced dimension", "space."))
  if (!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                  "partition", "pseudotime") | color_cells_by %in% 
                              names(colData(cds)), msg = paste("color_cells_by must one of", 
                                                               "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."))
    if (color_cells_by == "pseudotime") {
      tryCatch({
        pseudotime(cds, reduction_method = reduction_method)
      }, error = function(x) {
        stop(paste("No pseudotime for", reduction_method, 
                   "calculated. Please run order_cells with", 
                   "reduction_method =", reduction_method, "before attempting to color by pseudotime."))
      })
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers), 
                          msg = paste("Either color_cells_by or markers must", 
                                      "be NULL, cannot color by both!"))
  norm_method = match.arg(norm_method)
  group_cells_by = match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes), 
                          msg = paste("Either color_cells_by or genes must be", 
                                      "NULL, cannot color by both!"))
  if (show_trajectory_graph && is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize) {
    plotting_func <- ggrastr::geom_point_rast
  }
  else {
    plotting_func <- ggplot2::geom_point
  }
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    data_df$cell_group <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    stop("Error: unrecognized way of grouping cells.")
  }
  if (color_cells_by == "cluster") {
    data_df$cell_color <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    data_df$cell_color <- colData(cds)[data_df$sample_name, 
                                       color_cells_by]
  }
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = x, 
                                         prin_graph_dim_2 = y) %>% dplyr::mutate(sample_name = rownames(.), 
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[[reduction_method]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", 
                                                                     target = "to") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                           dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                          source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                         by = "source") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                                               dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                                              target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                                             by = "target")
  }
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2) {
      markers = unlist(genes[, 1], use.names = FALSE)
    }
    else {
      markers = genes
    }
    markers_rowData <- as.data.frame(subset(rowData(cds), 
                                            gene_short_name %in% markers | row.names(rowData(cds)) %in% 
                                              markers))
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), 
                                                     , drop = FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
      if (!is.null(dim(genes)) && dim(genes) >= 2) {
        genes = as.data.frame(genes)
        row.names(genes) = genes[, 1]
        genes = genes[row.names(cds_exprs), ]
        agg_mat = as.matrix(aggregate_gene_expression(cds, 
                                                      genes, norm_method = norm_method, scale_agg_values = FALSE))
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id", 
                                          "cell_id")
        if (is.factor(genes[, 2])) 
          markers_exprs$feature_id = factor(markers_exprs$feature_id, 
                                            levels = levels(genes[, 2]))
        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      }
      else {
        cds_exprs@x = round(10000 * cds_exprs@x)/10000
        markers_exprs = matrix(cds_exprs, nrow = nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id", 
                                          "cell_id")
        markers_exprs <- merge(markers_exprs, markers_rowData, 
                               by.x = "feature_id", by.y = "row.names")
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <- as.character(markers_exprs$feature_id)
        }
        else {
          markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        }
        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | 
                                                !as.character(markers_exprs$feature_label) %in% 
                                                markers, as.character(markers_exprs$feature_id), 
                                              as.character(markers_exprs$feature_label))
        markers_exprs$feature_label <- factor(markers_exprs$feature_label, 
                                              levels = markers)
        if (norm_method == "size_only") 
          expression_legend_label = "Expression"
        else expression_legend_label = "log10(Expression)"
      }
      if (scale_to_range) {
        markers_exprs = dplyr::group_by(markers_exprs, 
                                        feature_label) %>% dplyr::mutate(max_val_for_feature = max(value), 
                                                                         min_val_for_feature = min(value)) %>% dplyr::mutate(value = 100 * 
                                                                                                                               (value - min_val_for_feature)/(max_val_for_feature - 
                                                                                                                                                                min_val_for_feature))
        expression_legend_label = "% Max"
      }
    }
  }
  if (label_cell_groups && is.null(color_cells_by) == FALSE) {
    if (is.null(data_df$cell_color)) {
      if (is.null(genes)) {
        message(paste(color_cells_by, "not found in colData(cds), cells will", 
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }
    else {
      if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        if (label_groups_by_cluster && is.null(data_df$cell_group) == 
            FALSE) {
          text_df = data_df %>% dplyr::group_by(cell_group) %>% 
            dplyr::mutate(cells_in_cluster = dplyr::n()) %>% 
            dplyr::group_by(cell_color, add = TRUE) %>% 
            dplyr::mutate(per = dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                         text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>% 
            dplyr::top_n(labels_per_group, per)
        }
        else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>% 
            dplyr::mutate(per = 1)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                         text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>% 
            dplyr::top_n(labels_per_group, per)
        }
        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
      }
      else {
        message(paste("Cells aren't colored in a way that allows them to", 
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }
  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr, 
                                          value, NA))
    na_sub <- data_df[is.na(data_df$value), ]
    if (norm_method == "size_only") {
      g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
        plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                      stroke = I(cell_stroke), color = "grey95", 
                      alpha = alpha, data = na_sub) + plotting_func(aes(color = value), 
                                                                    size = I(cell_size), stroke = I(cell_stroke), 
                                                                    na.rm = TRUE) + viridis::scale_color_viridis(option = "viridis", 
                                                                                                                 name = expression_legend_label, na.value = "grey95", 
                                                                                                                 end = 0.8, alpha = alpha) + guides(alpha = FALSE) + 
        facet_wrap(~feature_label)
    }
    else {
      g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
        plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                      stroke = I(cell_stroke), color = "grey95", 
                      data = na_sub, alpha = alpha) + plotting_func(aes(color = log10(value + 
                                                                                        min_expr)), size = I(cell_size), stroke = I(cell_stroke), 
                                                                    na.rm = TRUE, alpha = alpha) + viridis::scale_color_viridis(option = "viridis", 
                                                                                                                                name = expression_legend_label, na.value = "grey95", 
                                                                                                                                end = 0.8, alpha = alpha) + guides(alpha = FALSE) + 
        facet_wrap(~feature_label)
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    if (color_cells_by %in% c("cluster", "partition")) {
      if (is.null(data_df$cell_color)) {
        g <- g + geom_point(color = I("gray"), size = I(cell_size), 
                            stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't", 
                      "color cells by cluster"))
      }
      else {
        g <- g + geom_point(aes(color = cell_color), 
                            size = I(cell_size), stroke = I(cell_stroke), 
                            na.rm = TRUE, alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by, 
                                           override.aes = list(size = 4)))
    }
    else if (class(data_df$cell_color) == "numeric") {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, 
                                            option = "C")
    }
    else {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by, 
                                           override.aes = list(size = 4)))
    }
  }
  if (show_trajectory_graph) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
                          color = I(trajectory_graph_color), linetype = "solid", 
                          na.rm = TRUE, data = edge_df)
    if (label_branch_points) {
      mst_branch_nodes <- monocle3:::branch_nodes(cds)
      branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes), 
                                                             sample_name)) %>% dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "white", fill = "black", size = I(graph_label_size * 
                                                                      1.5), na.rm = TRUE, branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                                  y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                                                                                                       size = I(graph_label_size), color = "white", 
                                                                                                                       na.rm = TRUE, branch_point_df)
    }
    if (label_leaves) {
      mst_leaf_nodes <- monocle3:::leaf_nodes(cds)
      leaf_df <- ica_space_df %>% dplyr::slice(match(names(mst_leaf_nodes), 
                                                     sample_name)) %>% dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "black", fill = "lightgray", size = I(graph_label_size * 
                                                                          1.5), na.rm = TRUE, leaf_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                              y = "prin_graph_dim_2", label = "leaf_idx"), 
                                                                                                                   size = I(graph_label_size), color = "black", 
                                                                                                                   na.rm = TRUE, leaf_df)
    }
    if (label_roots) {
      mst_root_nodes <- monocle3:::root_nodes(cds)
      root_df <- ica_space_df %>% dplyr::slice(match(names(mst_root_nodes), 
                                                     sample_name)) %>% dplyr::mutate(root_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "black", fill = "white", size = I(graph_label_size * 
                                                                      1.5), na.rm = TRUE, root_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                          y = "prin_graph_dim_2", label = "root_idx"), 
                                                                                                               size = I(graph_label_size), color = "black", 
                                                                                                               na.rm = TRUE, root_df)
    }
  }
  if (label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df, mapping = aes_string(x = "text_x", 
                                                                           y = "text_y", label = "label"), size = I(group_label_size))
    if (is.null(markers_exprs)) 
      g <- g + theme(legend.position = "none")
  }
  g <- g + monocle_theme_opts() + xlab(paste(reduction_method, 
                                             x)) + ylab(paste(reduction_method, y)) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  g
}


ad_plot_stacked_bar <- function(interval_vec, #vector of cell min_dist info
                                label_vec,    #vector of cell label(Top-level etc)
                                palette_list, #named list contains color and name of each cluster
                                title = "",
                                dilate_vec, img_size){
  
  tmp_df <- data.frame(label = rep(names(palette_list),each = length(dilate_vec)+1), 
                       interval = rep(c(names(dilate_vec),"overall"), length(palette_list)),
                       n = 0)
  tmp_df$label <- factor(tmp_df$label,levels = names(palette_list))
  tmp_df$n <- apply(tmp_df,1,
                    FUN = function(x){
                      if(x[["interval"]] == "overall") sum(label_vec == x[["label"]])
                      else sum(label_vec == x[["label"]] & interval_vec == x[["interval"]])
                    })
  for(i in names(dilate_vec)){
    tmp_df[tmp_df$interval == i,3] <- 10890000*tmp_df[tmp_df$interval == i,3]/dilate_vec[i]
  }
  tmp_df[tmp_df$interval == "overall",3] <- 10890000*tmp_df[tmp_df$interval == "overall",3]/img_size
  g <- ggplot(data = tmp_df,aes(x= interval, y=n, fill=label, order = label))+
    geom_bar(stat="identity") + 
    labs(title = title) + ylab("Density")+
    scale_fill_manual("",values = palette_list,
                      breaks = names(palette_list),
                      labels = names(palette_list)) 
  #theme_classic()
  
  return(g)
}

ad_s_gene_cl <- function(starmap_obj, cell_meta, sample_name, 
                         min_pct = 0.05, p_val = 0.01, seq_num = 6, p_filter = F){
  #Get Gene Expr in 0-40um filtered by pct
  #min_pct <- 0.05
  #starmap_obj <- starmap_raw
  #starmap_obj@misc <- list()
  #cell_meta <- starmap_raw@misc$ADmouse_9723_2_plaque_all$ADmouse_9723_2_cell_meta
  #ample_name <- "ADmouse_9723_2"
  #cell_meta <- starmap@misc$AD_mouse9494_plaque_all$AD_mouse9494_cell_meta
  #sample_name <- "AD_mouse9494"
  #p_val <- 0.01
  
  gene_list <- row.names(starmap_obj@assays[["RNA"]]@counts)[
    rowSums(starmap_obj@assays[["RNA"]]@counts[,starmap_obj@meta.data[["sample"]] == sample_name] > 0)
    > sum(starmap_obj@meta.data[["sample"]] == sample_name) * min_pct]
  
  tmp_df <- matrix(0, nrow = length(gene_list), ncol = 5) %>% as.data.frame(row.names = gene_list)
  colnames(tmp_df)[1:4] <- paste("R",1:4,"0",sep = "")
  colnames(tmp_df)[5] <- paste("R","40+",sep = "")
  #Raw Counts
  #tmp_mat <- starmap_obj@assays[["RNA"]]@counts[gene_list,starmap_obj@meta.data[["sample"]] == sample_name]
  #Scaled
  tmp_mat <- starmap_obj@assays[["scaled"]]@data[gene_list,starmap_obj@meta.data[["sample"]] == sample_name]
  for (i in 1:4) {
    tmp_df[,i] <- rowSums(tmp_mat[,cell_meta$interval == paste(i,"0um",sep = "")]) /
      sum(cell_meta$interval == paste(i,"0um",sep = ""))
  }
  tmp_df[,5] <- rowSums(tmp_mat[,cell_meta$interval %in% c("other","50um")]) / sum(cell_meta$interval %in% c("other","50um"))
  
  
  #kmeans try
  library(factoextra)
  fviz_nbclust(t(scale(t(as.matrix(tmp_df)))),kmeans,k.max = 15,method = "wss")
  set.seed(2021)
  tmp_df.km_res <- kmeans(t(scale(t(as.matrix(tmp_df)))), seq_num, nstart = 25)
  max_int_vec <- rep(0,seq_num)
  for(i in 1:seq_num){
    tmp_vec <- tmp_df[row.names(tmp_df) %in% gene_list[tmp_df.km_res$cluster == i],] %>% colSums() %>% .[1:4]
    max_int_vec[i] <- which(tmp_vec == max(tmp_vec))
  }
  names(max_int_vec) <- as.character(1:seq_num)
  #Stat test
  gene_df <- data.frame(gene = gene_list, cluster = tmp_df.km_res$cluster, max_int = 0, p_val = 1)
  gene_df$max_int <- apply(gene_df,1,FUN = function(x){max_int_vec[x[["cluster"]] ]})
  gene_df$p_val <- apply(gene_df,1,
                         FUN = function(x){
                           i = x[["max_int"]]
                           tmp_vec1 <- tmp_mat[x[["gene"]],cell_meta$interval == paste(i,"0um",sep = "")]
                           tmp_vec2 <- tmp_mat[x[["gene"]],cell_meta$interval != paste(i,"0um",sep = "")]
                           a <- wilcox.test(tmp_vec1,tmp_vec2, alternative = "greater")$p.value
                           return(a)
                         })
  
  v_heatmap <- NULL
  #row.names(tmp_df) <- str_to_title(row.names(tmp_df))
  if(p_filter == F){
    for(i in 1:seq_num){
      v_heatmap <- v_heatmap %v% Heatmap(t(scale(t(tmp_df[tmp_df.km_res$cluster == i, 1:5]))),
                                         cluster_columns = F,cluster_rows  = T, show_row_names = F, show_row_dend = F,
                                         row_title = paste("C",i),border_gp = gpar(col = "black", lty = 1, lwd = 2),
                                         column_title = paste(sample_name,"PCT:",min_pct,"P-val",p_val)) 
    }
  }else{
    for(i in 1:seq_num){
      v_heatmap <- v_heatmap %v% Heatmap(t(scale(t(tmp_df[tmp_df.km_res$cluster == i & 
                                                            gene_df$p_val < p_val, 1:5]))),
                                         cluster_columns = F,cluster_rows  = T, show_row_names = F, show_row_dend = F,
                                         row_title = paste("C",i),border_gp = gpar(col = "black", lty = 1, lwd = 2),
                                         column_title = paste(sample_name,"PCT:",min_pct,"P-val",p_val)) 
    }
  }
  
  #Save res
  gene_cl_res <- list()
  gene_cl_res[["max_int_vec"]] <- max_int_vec
  gene_cl_res[["p_val"]] <- p_val
  gene_cl_res[["expr_df"]] <- tmp_df
  gene_cl_res[["km_res"]] <- tmp_df.km_res
  gene_cl_res[["stat_test"]] <- gene_df
  gene_cl_res[["combined_heatmap"]] <- v_heatmap
  return(gene_cl_res)
}

ad_area <- function(dims_array){dims_array[1]*dims_array[2]}

#####p-Tau#####

ad_tau_grid_analysis <- function(tau_img,cell_meta,sample_name, 
                                 block_size, palette_list, level = "top_level", 
                                 plaque_img = NULL, region_img = NULL){
  plot_list <- list()
  tmp_img <- Image(tau_img!=0,colormode = Grayscale)
  tmp_mat <- matrix(0,nrow = floor(dim(tmp_img)[1] / block_size),
                    ncol = floor(dim(tmp_img)[2] / block_size))
  for(i in 1:dim(tmp_mat)[1]){
    for(j in 1:dim(tmp_mat)[2]){
      tmp_mat[i,j] <- sum(tmp_img[((i-1)*block_size+1):(i*block_size),((j-1)*block_size+1):(j*block_size)])
    }
  }
  require(ComplexHeatmap)
  require(circlize)
  plot_list[["tau_heatmap"]] <- Heatmap(tmp_mat,cluster_rows = F,cluster_columns = F,
                                        col =  colorRamp2(c(0, quantile(tmp_mat,probs = c(0.99))), 
                                                          c("white", "darkred")))
  #Cell Idty Stat
  cell_stat <- data.frame(x = rep(c(1:dim(tmp_mat)[1]), dim(tmp_mat)[2]),
                          y = rep(c(1:dim(tmp_mat)[2]), each = dim(tmp_mat)[1]))
  cell_stat[["tau_level"]] <- apply(cell_stat,1,FUN = function(x){tmp_mat[x[["x"]], x[["y"]] ]})
  for(i in names(palette_list)){
    cell_stat[[i]] <- apply(cell_stat,1,
                            FUN = function(x){
                              sum(cell_meta$x >= (x[["x"]]-1)*block_size+1 &
                                    cell_meta$x <= x[["x"]]*block_size &
                                    cell_meta$y >= (x[["y"]]-1)*block_size+1 &
                                    cell_meta$y <= x[["y"]]*block_size &
                                    cell_meta[[level]] == i)
                            })
    #print(i)
  }
  if(!is.null(region_img)){
    cell_stat[["region"]] <- apply(cell_stat,1,
                                   FUN = function(x){
                                     region_img[round((x[["x"]]-1/2)*block_size), round((x[["y"]]-1/2)*block_size)]
                                   })
  }
  require(reshape2)
  if(!is.null(plaque_img)){
    cell_stat[["plaque"]] <- apply(cell_stat,1,
                                   FUN = function(x){
                                     return(sum(plaque_img[((x[["x"]]-1)*block_size+1):(x[["x"]]*block_size),
                                                           ((x[["y"]]-1)*block_size+1):(x[["y"]]*block_size)]))
                                   })
    tmp_df <- data.frame(tau_interval = c("0%","50%","100%","100%","100%"),
                         min_tau = c(0,quantile(cell_stat$tau_level[cell_stat$tau_level != 0])[c(1,3,3,3)]),
                         max_tau = c(quantile(cell_stat$tau_level[cell_stat$tau_level != 0])[c(1,3)],rep(block_size^2,3)),
                         plaque = c("any","any","any","with_plaque","without_plaque"))
    for(i in 1:5){
      for(j in names(palette_list)){
        if(tmp_df$plaque[i] == "any"){
          tmp_df[[j]][i] <- sum(cell_stat[[j]][cell_stat$tau_level >= tmp_df$min_tau[i] &
                                                 cell_stat$tau_level < tmp_df$max_tau[i]])/
            sum(cell_stat$tau_level >= tmp_df$min_tau[i] & cell_stat$tau_level < tmp_df$max_tau[i])
        }else if(tmp_df$plaque[i] == "with_plaque"){
          tmp_df[[j]][i] <- sum(cell_stat[[j]][cell_stat$tau_level >= tmp_df$min_tau[i] &
                                                 cell_stat$tau_level < tmp_df$max_tau[i] & cell_stat$plaque > 0])/
            sum(cell_stat$tau_level >= tmp_df$min_tau[i] & cell_stat$tau_level < tmp_df$max_tau[i] & cell_stat$plaque > 0)
        }else if(tmp_df$plaque[i] == "without_plaque"){
          tmp_df[[j]][i] <- sum(cell_stat[[j]][cell_stat$tau_level >= tmp_df$min_tau[i] &
                                                 cell_stat$tau_level < tmp_df$max_tau[i] & cell_stat$plaque == 0])/
            sum(cell_stat$tau_level >= tmp_df$min_tau[i] & cell_stat$tau_level < tmp_df$max_tau[i] & cell_stat$plaque == 0)
        }
      }
    }
    tmp_df_bk <- tmp_df 
    tmp_df <- tmp_df %>% melt(id.vars = c("tau_interval","min_tau","max_tau","plaque"))
    tmp_df$tau_interval <- factor(tmp_df$tau_interval, levels = c("0%","50%","100%"))
    tmp_df$variable <- factor(tmp_df$variable, levels = names(palette_list))
    tmp_df[["title"]] <- factor(paste(tmp_df$tau_interval,tmp_df$plaque,sep = "_"),
                                levels = c("0%_any","50%_any","100%_any","100%_with_plaque","100%_without_plaque"))
    
    plot_list[["stacked_barplot"]] <-
      ggplot(data = tmp_df,aes(x = title, y = value, fill = variable, order = variable))+
      geom_bar(stat="identity") + 
      ggtitle(paste(sample_name,"Cell Density in Tau Free blocks and 0-25%/.../75-100% quantile blocks with tau")) +
      scale_fill_manual("",values = palette_list,
                        breaks = names(palette_list),
                        labels = names(palette_list)) +
      theme_classic()
    
  }else{
    #Ridge Barplot
    tmp_df <- data.frame(tau_interval = c("0%","50%","100%"),
                         #min_tau = c(0,1,10,29,68),
                         #max_tau = c(1,10,29,68,999),
                         min_tau = c(0,quantile(cell_stat$tau_level[cell_stat$tau_level != 0])[c(1,3)]),
                         max_tau = c(quantile(cell_stat$tau_level[cell_stat$tau_level != 0])[c(1,3)],block_size^2))
    
    for(i in 1:3){
      for(j in names(palette_list)){
        tmp_df[[j]][i] <- sum(cell_stat[[j]][cell_stat$tau_level >= tmp_df$min_tau[i] &
                                               cell_stat$tau_level < tmp_df$max_tau[i]])/
          sum(cell_stat$tau_level >= tmp_df$min_tau[i] & cell_stat$tau_level < tmp_df$max_tau[i])
      }
      
    }
    tmp_df_bk <- tmp_df 
    tmp_df <- tmp_df %>% melt(id.vars = c("tau_interval","min_tau","max_tau"))
    tmp_df$tau_interval <- factor(tmp_df$tau_interval, levels = c("0%","50%","100%"))
    tmp_df$variable <- factor(tmp_df$variable, levels = names(palette_list))
    
    plot_list[["stacked_barplot"]] <-
      ggplot(data = tmp_df,aes(x = tau_interval, y = value, fill = variable, order = variable))+
      geom_bar(stat="identity") + 
      ggtitle(paste(sample_name,"Cell Density in Tau Free blocks and 0-25%/.../75-100% quantile blocks with tau")) +
      scale_fill_manual("",values = palette_list,
                        breaks = names(palette_list),
                        labels = names(palette_list)) +
      theme_classic()
  }
  cell_stat[["tau_group"]] <- apply(cell_stat,1,
                                    FUN = function(x){
                                      tau_level <- as.numeric(x[["tau_level"]])
                                      plaque_level <- as.numeric(x[["plaque"]])
                                      
                                      tmp_vec <- quantile(cell_stat$tau_level[cell_stat$tau_level != 0], 
                                                          probs = c(0,0.5,1))
                                      if(tau_level < tmp_vec[1]) "Zero"
                                      else if(tau_level < tmp_vec[2]) "Low"
                                      else if(tau_level >= tmp_vec[2] & plaque_level == 0) "High w/o ABeta"
                                      else if(tau_level >= tmp_vec[2] & plaque_level > 0) "High w ABeta"
                                      else "NA"
                                    })
  
  
  return(list("plot_list" = plot_list,
              "cell_stat" = cell_stat,
              "interval_stat" = tmp_df,
              "interval_stat_wide" = tmp_df_bk))
}

ad_tau_stat <- function(cell_stat, cell_type_vec){
  tmp_df <- data.frame(cell_type = cell_type_vec)
  for(i in c("Low","High","High w ABeta","High w/o ABeta")){
    tmp_df[[i]] <- apply(tmp_df,1,
                         FUN = function(x){
                           if(i == "High") tmp_vec <- cell_stat[[x[["cell_type"]]]][cell_stat$tau_group %in% c("High w ABeta", "High w/o ABeta")]
                           else tmp_vec <- cell_stat[[x[["cell_type"]]]][cell_stat$tau_group == i]
                           tmp_obj <- t.test(x = tmp_vec,
                                             y = cell_stat[[x[["cell_type"]]]][cell_stat$tau_group == "Zero"],
                                             alternative = "greater")
                           tmp_obj$p.value
                         })
  }
  tmp_df
}

ad_combined_tau_density <- function(palette_list,stat_1,stat_2,title_text = ""){
  tmp_df <- data.frame(cell_type = stat_1$variable,
                       title = stat_1$title,
                       sample1 = stat_1$value,
                       sample2 = stat_2$value)
  tmp_df$average <- (tmp_df$sample1+tmp_df$sample2)/2
  tmp_df <- tmp_df %>% subset(subset = cell_type != "DG")
  palette_list <- palette_list[names(palette_list) != "DG"]
  tmp_df$title <- factor(tmp_df$title,levels = c("0%_any","50%_any","100%_any","100%_with_plaque","100%_without_plaque"))
  p <-
    ggplot(data = tmp_df,aes(x = title, y = average, fill = cell_type, order = cell_type))+
    geom_bar(stat="identity") + 
    ggtitle(paste("p-Tau Cell Density",title_text)) +
    scale_fill_manual("",values = palette_list ,
                      breaks = names(palette_list),
                      labels = names(palette_list)) +
    theme_classic()
}
