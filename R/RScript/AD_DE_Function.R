ad_D_vs_C_DE <- function(starmap_obj,
                         top_level_list,
                         time_name,
                         p_cutoff = 0.05,
                         logfc_cutoff = 0.1){
  
  tmp <- subset(starmap_obj,subset = top_level %in% top_level_list & time %in% time_name)
  tmp <- SetIdent(object = tmp, value  = tmp@meta.data[["group"]])
  require(tibble)
  tmp_df <- as.data.frame(AverageExpression(tmp)[[1]]) %>% rownames_to_column()
  marker_list <- FindMarkers(tmp,ident.1 = "disease",ident.2 = "control",logfc.threshold = logfc_cutoff) %>% rownames_to_column()
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
  marker_list <- FindMarkers(tmp,ident.1 = "disease",ident.2 = "control",logfc.threshold = 0) %>% rownames_to_column()
  
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
                              title = paste(top_level_list[1],time_name,"Disease vs. Control"),
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

#####GSEA#####
GSEA = function(gene_list, myGO, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600,
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  #print(dim(fgRes))
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
  
  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
  
  #print(dim(rbind(ups,downs)))
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ### Collapse redundant pathways
  Up = fgsea::collapsePathways(keepups, pathways = myGO, stats = gene_list,  nperm = 500, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns, myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  
  fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                              c( Up$mainPathways, Down$mainPathways))), ] %>% 
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                                 "Up-regulated" = "firebrick") ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  
  output = list("Results" = fgRes, "Plot" = g)
  return(output)
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



