#####Function Definition: Pseudo-time Analysis#####
psdtime_analysis <- function(starmap_obj, embed_df = NULL,
                             cell_type,
                             level = "top_level",
                             test = F,
                             umap_min_dist = 1,
                             euclidean_dist_ratio = 0.5,
                             geodesic_dist_ratio = 2/3,
                             min_branch_len = 10, npc = 10, nn = 20, ncenter = NULL, eps = NULL
){
  #level = "top_level"
  #cell_type = c("OPC","Oligo")
  #umap_min_dist = 1
  #embed_df <- embed_ls$Oligo_diff[,2:3]
  #For Micro Only
  #embed_df[,1] <- embed_df[,1]*100
  #embed_df[,2] <- embed_df[,2]*100
  #euclidean_dist_ratio = 0.5
  #geodesic_dist_ratio = 2/3
  #min_branch_len = 10
  tmp_seurat_obj <- starmap_obj
  tmp_seurat_obj@misc <- list()
  if(level == "top_level"){
    tmp_seurat_obj <- subset(tmp_seurat_obj, subset = top_level %in% c(cell_type))
  }
  expr_mat <- tmp_seurat_obj@assays[["RNA"]]@data
  cell_meta <- tmp_seurat_obj@meta.data
  colnames(expr_mat) <- paste(cell_meta$sample,cell_meta$orig_index,sep = "_")
  row.names(cell_meta) <- colnames(expr_mat)
  
  require(monocle3)
  cds <- new_cell_data_set(expr_mat,
                           cell_metadata = cell_meta,
                           gene_metadata = data.frame(gene_short_name =  row.names(tmp_seurat_obj@assays[["RNA"]]),
                                                      row.names = row.names(tmp_seurat_obj@assays[["RNA"]]),
                                                      highly_varible = tmp_seurat_obj@assays[["RNA"]]@meta.features[["highly_variable"]]))
  print("Preprocessing...")
  cds <- preprocess_cds(cds, num_dim = npc, norm_method = "none", 
                        use_genes = row.names(tmp_seurat_obj@assays[["RNA"]])[tmp_seurat_obj@assays[["RNA"]]@meta.features[["highly_variable"]] == T],
                        scaling = F)
  #Grid Parameter adjusting   
  #OPC/Oligo 1.25  Astro 1.25 Micro 0.50 
  #Ex        1.25  Inhi  0.50 CA1-3 0.50
  cds <- reduce_dimension(cds,max_components = 2,
                          umap.min_dist = umap_min_dist, umap.n_neighbors = nn)
  print("Clustering...")
  if(!is.null(embed_df)){
    #To ELIMINATE THE COLNAMES OF embed_df
    embed_df <- matrix(unlist(embed_df), nrow = dim(embed_df)[1], ncol = dim(embed_df)[2])
    row.names(embed_df) <- row.names(cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]])
    cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- embed_df
  }
  cds <- cluster_cells(cds)
  
  #OPC/Oligo  0.5 2/3 6
  #Astro      0.5 2/3 15
  #Micro      0.5 2/3 15
  #Ex         0.5 2/3 10
  #Inhi       0.5 2/3 10
  #CA1-3      0.5 2/3 10
  cds <- learn_graph(cds,close_loop = F,
                     learn_graph_control = list(euclidean_distance_ratio = euclidean_dist_ratio,
                                                geodesic_distance_ratio = geodesic_dist_ratio,
                                                minimal_branch_len = min_branch_len,
                                                orthogonal_proj_tip = F, ncenter = ncenter, eps = eps))
  
  cds <- order_cells(cds)
  plot_list <- list()
  plot_list[["pseudotime"]] <- plot_cells(cds,
                                          color_cells_by = "pseudotime",
                                          label_cell_groups=FALSE,
                                          label_leaves=FALSE,
                                          label_branch_points=FALSE,
                                          graph_label_size=1.5,
                                          cell_size = 1.2)
  
  plot_list[["cell_type"]] <- plot_cells(cds,
                                         color_cells_by = "cell_type",
                                         label_groups_by_cluster=FALSE,
                                         label_leaves=FALSE,
                                         label_branch_points=FALSE,
                                         group_label_size=8,
                                         cell_size = 1.2)+ theme(legend.position = "right")
  
  plot_list[["sample"]] <- plot_cells(cds,
                                      color_cells_by = "sample",
                                      label_groups_by_cluster=FALSE,
                                      label_leaves=FALSE,
                                      label_branch_points=FALSE,
                                      group_label_size=0,
                                      cell_size = 1.2)+ theme(legend.position = "right")
  
  tmp_df <- data.frame(tmp_seurat_obj@meta.data,
                       pseudotime = cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]],
                       UMAP1 = cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,1],
                       UMAP2 = cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]][,2])
  return_list <- list("Plot_list" = plot_list,
                      "cell_meta" = tmp_df,
                      "cds" = cds)
}


#Plot psdtime vs distance
ad_plot_psdtime_distr <- function(cell_meta, psdtime_res, sample, cell_type){
  tmp_df <- cell_meta %>% 
    .[.$top_level %in% cell_type,] %>%
    data.frame(#dist_group = pmin(floor(.$min_border_dist/33),5),
               pseudotime = psdtime_res[psdtime_res$sample %in% sample & psdtime_res$top_level %in% cell_type,c("pseudotime")])
  tmp_df$pseudotime <- tmp_df$pseudotime/max(tmp_df$pseudotime)
  tmp_df <- tmp_df %>% mutate(interval = "mean") %>% rbind(subset(tmp_df, interval != "other"))
  #tmp_df <- tmp_df %>% mutate(dist_group = "mean") %>% 
  #  rbind(subset(tmp_df,dist_group < 5)) %>%
  #  mutate(dist_group = as.factor(dist_group))
  #levels(tmp_df$dist_group) <- c(paste(1:5,"0um",sep = ""),"mean")
  
  tmp_df %>% ggplot( aes(x=interval, y=pseudotime, fill=interval)) +
    geom_boxplot() +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_bw() +
    theme(
      legend.position="none",
      panel.background = element_blank(),
      plot.title = element_text(size=11)
    ) + 
    scale_fill_manual(values= colorRamp(c("#7703fc", "#8403fc", "#be03fc", "#fc6b03", "#fceb03"),space = "rgb")#Based on mean value of each group to assign a color
                      #(aggregate(.~ dist_group, tmp_df[,13:14], median)$pseudotime) %>% `/`(255) %>% rgb()) +
                      (aggregate(.~ interval, tmp_df[,13:14], median)$pseudotime) %>% `/`(255) %>% rgb()) +
    ggtitle(paste(paste(sample),"Pseudotime Distribution",cell_type)) +
    xlab("")
}

#Highlight cells with graph
ad_psdtime_highlight <- function(meta_df, highlight_col, cds = NULL, show_trajectory_graph = F, title){
  plot_list <- list()
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>% 
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = 1, 
                                         prin_graph_dim_2 = 2) %>% dplyr::mutate(sample_name = rownames(.), 
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[["UMAP"]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
      dplyr::select_(source = "from",target = "to") %>% 
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                       by = "source") %>% 
      dplyr::left_join(ica_space_df %>% dplyr::select_(target = "sample_name", 
                                                       target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                       target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                       by = "target")
  }
  meta_df[["pseudotime"]] <- meta_df[["pseudotime"]]/max(meta_df[["pseudotime"]])
  tmp_col <- cell_type_palette_list[[title]]
  names(tmp_col) <- cell_type_list[[title]]
  for(i in c("all",levels(droplevels(meta_df[[highlight_col]])))){
    g <- ggplot() +
      theme(line = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      labs(title = paste(i,title))
    if(i != "all"){
      g <- g + geom_point(data = meta_df[meta_df[[highlight_col]] != i,], 
                          aes(x = UMAP1, y = UMAP2), colour = "gray90",alpha = 0.8)
    }
    if(highlight_col == "sample"){
      if(i == "all"){
        g <- g + geom_point(data = meta_df, shape = 16, 
                            aes(x = UMAP1, y = UMAP2, colour = pseudotime),alpha = 0.4, size = 3)
      }else{
        g <- g + geom_point(data = meta_df[meta_df[[highlight_col]] == i,], shape = 16, 
                            aes(x = UMAP1, y = UMAP2, colour = pseudotime),alpha = 0.4, size = 3)
      }
      g <- g + scale_color_gradientn(colours = c("#7703fc", "#8403fc", "#be03fc", "#fc6b03", "#fceb03"),
                                     values = c(0, 0.2, 0.4, 0.6, 0.8, 1))
      
    }else if(highlight_col == "cell_type"){
      
      if(i == "all"){
        g <- g + geom_point(data = meta_df, shape = 16, 
                            aes(x = UMAP1, y = UMAP2, color = cell_type),alpha = 0.4, size = 3)
      }else{
        g <- g + geom_point(data = meta_df[meta_df[[highlight_col]] == i,], shape = 16, 
                            aes(x = UMAP1, y = UMAP2, color = cell_type),alpha = 0.4, size = 3)
      }
      g <- g + scale_color_manual(values = tmp_col)
    }
    # colour = cell_type_palette_list[[title]][[which(cell_type_list[[title]] == i)]]) +#colour = pseudotime
    #scale_color_gradientn(colours = c("#7703fc", "#8403fc", "#be03fc", "#fc6b03", "#fceb03"),
    #                      values = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    
    #labs(title = paste(i,title,"Pseudotime"))
    if(show_trajectory_graph){
      mst_root_nodes <- cds@principal_graph_aux@listData[["UMAP"]][["root_pr_nodes"]]
      root_df <- ica_space_df[mst_root_nodes,]
      
      g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                       y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                       yend = "target_prin_graph_dim_2"), size = 1.0, 
                            color = "gray28", linetype = "solid", 
                            na.rm = TRUE, data = edge_df) + 
        geom_point(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = 0.75, 
                   color = "black", fill = "white", size = 2 * 1.5, na.rm = TRUE, root_df) + 
        geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "1"), size = 2, 
                  color = "black", na.rm = TRUE, root_df)
    }
    plot_list[[i]] <- g
  }
  
  return(plot_list)
}

ad_psdtime_gene_highlight <- function(meta_df,cds,gene_list,show_trajectory_graph = T){
  plot_list <- list()
  meta_df[["pseudotime"]] <- meta_df[["pseudotime"]]/max(meta_df[["pseudotime"]])
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>% 
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = 1, 
                                         prin_graph_dim_2 = 2) %>% dplyr::mutate(sample_name = rownames(.), 
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[["UMAP"]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
      dplyr::select_(source = "from",target = "to") %>% 
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                       by = "source") %>% 
      dplyr::left_join(ica_space_df %>% dplyr::select_(target = "sample_name", 
                                                       target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                       target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                       by = "target")
  }
  for(i in gene_list){
    plot_list[[i]] <- ggplot() +
      theme(line = element_blank(),
            panel.background = element_blank(),
            legend.position = "none") +
      geom_point(data = meta_df[cds@assays@data@listData[["counts"]][i,] <= 0,], 
                 aes(x = UMAP1, y = UMAP2), colour = "gray90",alpha = 0.8) +
      geom_point(data = meta_df[cds@assays@data@listData[["counts"]][i,] > 0,], shape = 16, 
                 aes(x = UMAP1, y = UMAP2,  colour = pseudotime),alpha = 0.4, size = 3) +
      scale_color_gradientn(colours = c("#7703fc", "#8403fc", "#be03fc", "#fc6b03", "#fceb03"),
                            values = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      
      labs(title = paste(i,"Pseudotime"))
    if(show_trajectory_graph){
      mst_root_nodes <- cds@principal_graph_aux@listData[["UMAP"]][["root_pr_nodes"]]
      root_df <- ica_space_df[mst_root_nodes,]
      
      plot_list[[i]] <- plot_list[[i]] + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                                                 y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                                                 yend = "target_prin_graph_dim_2"), size = 1.0, 
                                                      color = "gray28", linetype = "solid", 
                                                      na.rm = TRUE, data = edge_df) + 
        geom_point(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = 0.75, 
                   color = "black", fill = "white", size = 2 * 1.5, na.rm = TRUE, root_df) + 
        geom_text(aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "1"), size = 2, 
                  color = "black", na.rm = TRUE, root_df)
    }
  }
  return(plot_list)
}

