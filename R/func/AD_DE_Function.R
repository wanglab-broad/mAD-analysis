#####Function Definition: Differential Expression Analysis Function#####
ad_D_vs_C_DE <- function(starmap_obj,
                         cell_type_list,
                         time_name,
                         p_cutoff = 0.05,
                         logfc_cutoff = 0.1, sub_type = F){
  if(sub_type){
    tmp <- subset(starmap_obj,subset = cell_type %in% cell_type_list & time %in% time_name)
  }else{
    tmp <- subset(starmap_obj,subset = top_level %in% cell_type_list & time %in% time_name)
  }
 
  tmp <- SetIdent(object = tmp, value  = tmp@meta.data[["group"]])
  require(tibble)
  tmp_df <- as.data.frame(AverageExpression(tmp)[[1]]) %>% rownames_to_column()
  marker_list <- FindMarkers(tmp,ident.1 = "disease",ident.2 = "control",logfc.threshold = logfc_cutoff,min.pct = 0.05) %>% rownames_to_column()
  marker_list <- marker_list[marker_list$p_val < p_cutoff ,]
  #write.csv(marker_list,paste("disease_vs_ctrl",top_level_list[i],".csv",sep = "_"),row.names = T)
  return(marker_list)
}

ad_D_vs_C_DE_volcano <- function(starmap_obj,
                                 top_level_list,
                                 time_name,
                                 fig_title = NULL,
                                 FC_cutoff = 0.1,
                                 p_cutoff = 0.05,
                                 Pcap = NULL,
                                 FCcap = NULL,
                                 draw_connector = F,
                                 xrange = NULL,
                                 yrange = NULL){
  
  tmp <- subset(starmap_obj,subset = top_level %in% top_level_list & time %in% time_name)
  tmp <- SetIdent(object = tmp, value  = tmp@meta.data[["group"]])
  require(tibble)
  tmp_df <- as.data.frame(AverageExpression(tmp)[[1]]) %>% rownames_to_column()
  marker_list <- FindMarkers(tmp,ident.1 = "disease",ident.2 = "control",logfc.threshold = 0, min.pct = 0.05) %>% rownames_to_column()
  
  tmp_df <- left_join(marker_list,tmp_df,by = "rowname")
  require(EnhancedVolcano)
  if(!is.null(Pcap)) marker_list$p_val[marker_list$p_val < Pcap] <- Pcap
  if(!is.null(FCcap)){
    marker_list$avg_logFC[marker_list$avg_logFC > FCcap] <- FCcap
    marker_list$avg_logFC[marker_list$avg_logFC < -FCcap] <- -FCcap
  }
  tmp_col <- ifelse(marker_list$p_val>p_cutoff,"grey",
                    ifelse(abs(marker_list$avg_logFC) < FC_cutoff,
                           "darkgreen",ifelse(marker_list$avg_logFC < 0, "darkblue", "darkred")))
  require(stringr)
  row.names(marker_list) <- str_to_title(row.names(marker_list))
  names(tmp_col) <- row.names(marker_list)
  marker_list$rowname <- str_to_title(marker_list$rowname)
  tmp_plot <- EnhancedVolcano(marker_list,
                              lab = marker_list$rowname,
                              x="avg_logFC",y="p_val",
                              title = paste(fig_title,time_name,"Disease vs. Control"),
                              subtitle = " ",
                              xlim = ifelse(c(is.null(xrange),is.null(xrange)),
                                            c(-max(abs(marker_list$avg_logFC)) - 0.1,
                                              max(abs(marker_list$avg_logFC)) + 0.1), xrange),
                              ylim = ifelse(c(is.null(yrange),is.null(yrange)),
                                            c(0, max(-log10(marker_list$p_val), na.rm=TRUE) + 0.1)
                                            ,yrange),
                              FCcutoff = FC_cutoff,pCutoff = p_cutoff,
                              colCustom = tmp_col,
                              #Disable caption, grid and lines
                              caption = "",
                              captionLabSize = 0,
                              legendPosition = "none",
                              gridlines.major = F,
                              gridlines.minor = F,
                              drawConnectors = draw_connector,
                              widthConnectors = 0.5,
                              endsConnectors = 'none',
                              colConnectors = 'grey50',
                              arrowhead = F,
                              shape = 16,
                              pointSize = 1.5,
                              colAlpha = 0.9)
  
  return(tmp_plot)
}

#Median Scaling
#Divided by median
ad_starmap_normalize <- function(starmap_obj){
  starmap_scaled <- starmap_obj
  median_vec <- rep(0,length(levels(starmap_scaled@meta.data[["sample"]])))
  names(median_vec) <- levels(starmap_scaled@meta.data[["sample"]])
  for(i in levels(starmap_scaled@meta.data[["sample"]])){
    median_vec[[i]] <- median(starmap_scaled@meta.data[["n_counts"]][starmap_scaled@meta.data[["sample"]] == i])
  }
  
  scale_factor <- mean(na.omit(median_vec))
  
  for(i in levels(starmap_scaled@meta.data[["sample"]])){
    starmap_scaled@assays[["RNA"]]@data[,starmap_scaled@meta.data[["sample"]] == i] <- 
      log2(scale_factor*starmap_scaled@assays[["RNA"]]@data[,starmap_scaled@meta.data[["sample"]] == i]/
             median(starmap_scaled@meta.data[["n_counts"]][starmap_scaled@meta.data[["sample"]] == i])+1)
  }
  starmap_scaled@misc <- list()
  return(starmap_scaled)
}
#Median Scaling
#Divided by median for each cell type
ad_cell_type_median_norm <- function(starmap_obj){
  starmap_scaled <- starmap_obj
  median_vec <- rep(0,length(levels(starmap_scaled@meta.data[["sample"]])))
  names(median_vec) <- levels(starmap_scaled@meta.data[["sample"]])
  for(cell_type in levels(starmap_scaled@meta.data[["top_level"]])){
    for(i in levels(starmap_scaled@meta.data[["sample"]])){
      median_vec[[i]] <- median(starmap_scaled@meta.data[["n_counts"]][starmap_scaled@meta.data[["sample"]] == i & starmap_scaled@meta.data[["top_level"]] == cell_type])
    }
    
    scale_factor <- mean(na.omit(median_vec))
    
    for(i in levels(starmap_scaled@meta.data[["sample"]])){
      starmap_scaled@assays[["RNA"]]@data[,starmap_scaled@meta.data[["sample"]] == i  & starmap_scaled@meta.data[["top_level"]] == cell_type] <- 
        scale_factor*starmap_scaled@assays[["RNA"]]@data[,starmap_scaled@meta.data[["sample"]] == i & starmap_scaled@meta.data[["top_level"]] == cell_type]/
               median(starmap_scaled@meta.data[["n_counts"]][starmap_scaled@meta.data[["sample"]] == i & starmap_scaled@meta.data[["top_level"]] == cell_type])
    }
  }
  
  starmap_scaled@misc <- list()
  return(starmap_scaled)
}

