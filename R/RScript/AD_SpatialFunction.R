#####Function Definition: plot#####
ad_plot_cell_distr <- function(min_dist_vec, #vector of cell min_dist info
                               label_vec,    #vector of cell label(Top-level etc)
                               scale = "none", #none/percentage/density
                               palette_list, #named list contains color and name of each cluster
                               title = "",
                               step = 10, 
                               rings = 5, 
                               area_vec = NULL,     #length = rings, required when scale = "density
                               img_size = NULL, 
                               overall = T){
  for(i in 1:rings){
    if(i == 1){
      tmp_df <- data.frame(Dist_group = paste(i*step,"um",sep = ""),
                           as.data.frame(table(label_vec[min_dist_vec < 3.3*step & min_dist_vec >= 0])))
    }else{
      tmp_df <- rbind(tmp_df,
                      data.frame(Dist_group = paste(i*step,"um",sep = ""),
                                 as.data.frame(table(label_vec[min_dist_vec>= 3.3*step *(i-1) & min_dist_vec < 3.3*step*i]))))
    }
  }
  colnames(tmp_df) <- c("Dist_group","label","Freq")
  if(overall){
    tmp_df <- rbind(tmp_df,
                    data.frame(Dist_group = "overall",
                               as.data.frame(table(label_vec, dnn = list("label")), responseName = "Freq")))
  }
  
  tmp_df <- tmp_df[as.character(tmp_df$label) %in% names(palette_list),]
  if(scale == "percentage"){
    for(i in 1:rings){
      tmp_df[tmp_df$Dist_group == paste(i*step,"um",sep = ""),3] = 
        tmp_df[tmp_df$Dist_group == paste(i*step,"um",sep = ""),3]/sum(tmp_df[tmp_df$Dist_group == paste(i*step,"um",sep = ""),3])
    }
    if(overall){
      tmp_df[tmp_df$Dist_group == "overall",3] = 
        tmp_df[tmp_df$Dist_group == "overall",3]/sum(tmp_df[tmp_df$Dist_group == "overall",3])
    }
    
  }else if(scale == "density"){
    for(i in 1:rings){
      tmp_df[tmp_df$Dist_group == paste(i*step,"um",sep = ""),3] = 1000000 *
        tmp_df[tmp_df$Dist_group == paste(i*step,"um",sep = ""),3]/area_vec[i]*3.3^2
    }
    if(overall){
      tmp_df[tmp_df$Dist_group == "overall",3] = 1000000 *
        tmp_df[tmp_df$Dist_group == "overall",3]/img_size*3.3^2
    }
  }
  g <- ggplot(data = tmp_df,aes(x= Dist_group, y=Freq, fill=label, order = label))+
    geom_bar(stat="identity") + 
    labs(title = title) +
    scale_fill_manual("",values = palette_list,
                      breaks = names(palette_list),
                      labels = names(palette_list))
  if(overall){
    g <- g + scale_x_discrete(limits = c(paste(c(1:rings)*step,"um",sep = ""),"overall")) 
  }
  
  return(g)
}

ad_hl_cell_distr <- function(min_dist_vec, 
                             label_vec,   
                             scale = "none", #with other cell type colored in gray50
                             palette_list, 
                             title = "",
                             step = 10, 
                             rings = 5,
                             area_vec = NULL,     
                             img_size = NULL, 
                             overall = T,
                             y_min = 0, 
                             y_max = 1){
  
  palette_list <- c("Other" = "#AFAFAF", palette_list)
  label_vec <- as.character(label_vec)
  label_vec[!label_vec %in% names(palette_list)] <- "Other"
  label_vec <- factor(label_vec, levels = names(palette_list))
  
  g <- ad_plot_cell_distr(min_dist_vec,label_vec,scale,palette_list,title,step, rings, area_vec, img_size, overall)
  if(scale != "density"){
    g + coord_cartesian(ylim = c(y_min,y_max))
  }
  
}
#####Function Definition: preprocess and min_dist calc#####
#TEST
#starmap_obj       <- starmap
#img_list          <- starmap_obj@misc$AD_mouse9494_morph
#starmap_obj@misc  <- list()
#sample_name       <- "AD_mouse9494"
#min_area          <- 400
#region            <- "Cortex"
#palette_list <- as.character(unlist(starmap_raw@misc[["top_hex_dict"]]))[colnames(starmap_raw@misc[["top_hex_dict"]]) != "DG"]
#names(palette_list) <- colnames(starmap_raw@misc[["top_hex_dict"]])[colnames(starmap_raw@misc[["top_hex_dict"]]) != "DG"]
#rm(list = c("starmap_obj","img_list","cell_meta","plaque_meta","sample_name","min_area","region"))        
#
ad_calc_plaque_dist <- function(starmap_obj,
                                img_list,
                                min_area = 400,
                                sample_name,
                                palette_list,
                                region = "all"){
  #Calc Plaque Center
  print("Load Image")
  plaque <- Image(img_list[["plaque"]],colormode = Grayscale) %>% bwlabel()
  print("Compute Plaque Features")
  plaque_meta <- cbind(computeFeatures.moment(plaque),computeFeatures.shape(plaque)) %>% 
    as.data.frame()  %>% data.frame(., plaque_id = as.numeric(row.names(.)))
  
  cell_meta <- data.frame(starmap_obj@meta.data[starmap_obj@meta.data[["sample"]] == sample_name,],
                          min_center_dist = sqrt(sum(dim(plaque)^2)),
                          min_border_dist = sqrt(sum(dim(plaque)^2)),
                          nearest_plaque = 0)
  
  #Organize cell_meta
  #orig_index |x  |y  |batch  |time |group  |top_level  |cell_type  |region |min_dist |
  if(is.null(cell_meta$region)){
    cell_meta <-data.frame(orig_index = cell_meta$orig_index,
                           x = cell_meta$x, y = cell_meta$y,
                           batch = cell_meta$batch, time = cell_meta$time, group = cell_meta$group,
                           top_level = cell_meta$top_level, cell_type = cell_meta$cell_type,
                           region = NA,
                           min_center_dist = cell_meta$min_center_dist,
                           min_border_dist = cell_meta$min_border_dist,
                           nearest_plaque = 0)
  }else{
    cell_meta <-data.frame(orig_index = cell_meta$orig_index,
                           x = cell_meta$x, y = cell_meta$y,
                           batch = cell_meta$batch, time = cell_meta$time, group = cell_meta$group,
                           top_level = cell_meta$top_level, cell_type = cell_meta$cell_type,
                           region = cell_meta$region,
                           min_center_dist = cell_meta$min_center_dist,
                           min_border_dist = cell_meta$min_border_dist,
                           nearest_plaque = 0)
  }
  
  #Region Selection
  if(region != "all"){
    plaque_meta <- data.frame(plaque_meta,region = -1)
    plaque_meta$region <- apply(plaque_meta,1,
                                FUN = function(x){
                                  switch(img_list[["region"]][round(as.numeric(x[1])),round(as.numeric(x[2]))],
                                         '1' = "Cortex", '2' = "White Matter", '3' = "Hippocampus")
                                })
    if(sum(is.na(cell_meta$region)==F) == 0){
      cell_meta$region <- apply(cell_meta,1,
                                FUN = function(x){
                                  switch(img_list[["region"]][round(as.numeric(x[2])),round(as.numeric(x[3]))],
                                         '1' = "Cortex", '2' = "White Matter", '3' = "Hippocampus")
                                })
    }
  }
  #Save Raw palque meta for dilate
  plaque_meta_bk <- plaque_meta
  plaque_meta <- plaque_meta[plaque_meta$s.area > min_area,]
  
  if(region == "Sub-cortical"){
    cell_meta <- cell_meta[cell_meta$region %in% c("Hippocampus","White Matter"),]
    plaque_meta <- plaque_meta[plaque_meta$region %in% c("Hippocampus","White Matter"),]
  }else if(region != "all"){
    cell_meta <- cell_meta[cell_meta$region == region,]
    plaque_meta <- plaque_meta[plaque_meta$region == region,]
  }
  
  #Calc Min Dist
  print("Calc Min Plaque Dist")
  
  for(i in 1:dim(plaque_meta)[1]){
    tmp_list <- apply(cell_meta,MARGIN = 1,FUN = function(x){
      tmp_dist <- sqrt((as.numeric(x[["x"]]) - plaque_meta$m.cx[i])^2+(as.numeric(x[["y"]]) - plaque_meta$m.cy[i])^2)
      if(tmp_dist < as.numeric(x[["min_center_dist"]])){
        return(list(min_dist = tmp_dist, nearest_plaque = plaque_meta$plaque_id[i]))
      }else{
        return(list(min_dist = as.numeric(x[["min_center_dist"]]), nearest_plaque = as.numeric(x[["nearest_plaque"]])))
      }
    })
    #Update cell_meta
    tmp_list <- data.frame(matrix(unlist(tmp_list), nrow = length(tmp_list), byrow = T))
    cell_meta$min_center_dist <- as.numeric(tmp_list[,1])
    cell_meta$nearest_plaque <- as.numeric(tmp_list[,2])
  }
  
  #Get Border Distance
  print("Calc Min Dist to Plaque Border")
  cell_meta$min_border_dist <- apply(cell_meta,MARGIN = 1,
                                     function(x){
                                       cell_posX <- as.numeric(x[["x"]])
                                       cell_posY <- as.numeric(x[["y"]])
                                       min_center_dist <- as.numeric(x[["min_center_dist"]])
                                       nearest_plaque <- as.numeric(x[["nearest_plaque"]])
                                       j = which(plaque_meta$plaque_id == nearest_plaque)
                                       curr_dist <- max(0,min_center_dist - plaque_meta$s.radius.max[j] - 3) 
                                       #Move away from border
                                       ratioX <- (cell_posX - plaque_meta$m.cx[j])/min_center_dist
                                       ratioY <- (cell_posY - plaque_meta$m.cy[j])/min_center_dist
                                       
                                       posX <- round(cell_posX - curr_dist * ratioX)
                                       posY <- round(cell_posY - curr_dist * ratioY)
                                       cycle_count <- 0
                                       while(plaque[posX, posY] != nearest_plaque){
                                         curr_dist <- curr_dist + 1
                                         cycle_count <- cycle_count + 1
                                         posX <- round(cell_posX - curr_dist * ratioX)
                                         posY <- round(cell_posY - curr_dist * ratioY)
                                         if(cycle_count >= 1000){
                                           curr_dist <- -1
                                           break
                                         }
                                       }
                                       return(curr_dist)
                                     })
  
  #Calc Plaque 5-round dilate area 
  print("Calc size of dilated areas")
  if(region == "Sub-cortical"){
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area | 
                                                          !(plaque_meta_bk$region %in% c("Hippocampus","White Matter"))])
  }else if(region != "all"){
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area | 
                                                          plaque_meta_bk$region != region])
  }else{
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area])
  }
  
  plaque <- Image(plaque != 0,colormode = Grayscale)

  if(region == "all"){
    img_size <- dim(plaque)[1] * dim(plaque)[2]
  }else if(region == "Sub-cortical"){
    img_size <- sum(img_list[["region"]] == 2) + 
      sum(img_list[["region"]] == 3)
    tmp_img <- Image(img_list[["region"]] == 2 |img_list[["region"]] == 3)
  }else{
    x <- switch(region,
                "Cortex" = 1,"White Matter" = 2, "Hippocampus" = 3)
    img_size <- sum(img_list[["region"]] == x)
    tmp_img <- Image(img_list[["region"]] == x)
  }
  
  #Assign interval
  cell_meta <- data.frame(cell_meta,interval = "other")
  
  size_vec  <- rep(0,5); names(size_vec) <- paste(1:5,"0um",sep = "")
  dilate_img <- plaque
  dilate_img_bk <- plaque
  for(i in 1:5){
    dilate_img_bk <- dilate_img
    dilate_img <- dilate(plaque, kern = makeBrush(round(3.175*20*i), shape = "disc"))
    if(region != "all") dilate_img <- dilate_img*tmp_img
    if(i != 1) dilate_img2 <- dilate_img * (1 - dilate_img_bk)
    else dilate_img2 <- dilate_img
    size_vec[i] <- sum(dilate_img2) 
    
    cell_meta[["interval"]] <- apply(cell_meta,1,
                                     FUN = function(x){
                                       if(dilate_img2[as.numeric(x[["x"]]),as.numeric(x[["y"]])] != 0 ) return(paste(i,"0um",sep = ""))
                                       else return(x[["interval"]])
                                     })
    #print(table(cell_meta[["interval"]]))
  }
  
  #size_vec[1] <- tmp_vec[1] - sum(plaque)
  
  return_list <- list()
  return_list[[paste(sample_name,"plaque_dilate_area",sep = "_")]] <- size_vec
  return_list[[paste(sample_name,"cell_meta",sep = "_")]] <- cell_meta
  return_list[[paste(sample_name,"plaque_meta",sep = "_")]] <- plaque_meta
  return_list[[paste(sample_name,"plaque_img",sep = "_")]] <- plaque
  #Filter Palette list 
  #if(region == "Sub-cortical"){
  #  palette_list <- palette_list[names(palette_list) != "Ex"]
  #}else if(region == "Cortex"){
  #  palette_list <- palette_list[!names(palette_list) %in% c("CA1","CA2","CA3","DG")]
  #}else if(region == "Hippocampus"){
  #  palette_list <- palette_list[names(palette_list) != "Ex"] 
  #}else if(region == "White Matter"){
  #  palette_list <- palette_list[!names(palette_list) %in% c("CA1","CA2","CA3","DG","Ex","Inhi")]
  #}
  #plot_list <- list()
  #min_dist_vec, #vector of cell min_dist info
  #label_vec,    #vector of cell label(Top-level etc)
  #scale = "none", #none/percentage/density
  #palette_list, #named list contains color and name of each cluster
  #title = "",
  #step = 10, 
  #rings = 5, 
  #area_vec = NULL,     #length = rings, required when scale = "density
  #img_size = NULL, 
  #overall = T
  #plot_list[[1]] <- ad_plot_cell_distr(cell_meta$min_border_dist,cell_meta$top_level,scale = "none", palette_list,
  #                                     title = paste(sample_name,"Raw Counts"), overall = F)
  #plot_list[[2]] <- ad_plot_cell_distr(cell_meta$min_border_dist,cell_meta$top_level,scale = "density", palette_list,
  #                                     title = paste(sample_name,"Cells per sq. mm"),area_vec = size_vec,img_size = img_size, overall = T)
  #plot_list[[3]] <- ad_plot_cell_distr(cell_meta$min_border_dist,cell_meta$top_level,scale = "percentage",palette_list,
  #                                     title = paste(sample_name,"Percentage"))
  require(ggpubr)
  #tmp_plot <- ggarrange(plotlist = plot_list,nrow = 1,ncol = 1)
  #return_list[[paste(sample_name,"top_level_distr_plot",sep = "_")]] <- tmp_plot
  #return_list[[paste(sample_name,"t_l_d_plotlist",sep = "_")]] <- plot_list
  starmap_obj@misc[[paste(sample_name,"plaque",region,sep = "_")]] <- return_list
  return(starmap_obj)
}


#####Function Definition: DE Analysis#####
#####SpatialDE 25um Micro#####
ad_spatial_DE <- function(starmap_obj,
                          sample_name,
                          region = "all",
                          top_level = "all",
                          dist_cutoff = 25,
                          p_cutoff = 0.05){
  
  tmp <- subset(starmap_obj,subset = sample == sample_name)
  tmp <- AddMetaData(tmp, metadata = data.frame(min_border_dist = starmap_obj@misc[[paste(sample_name,"plaque",region,sep = "_")]]
                                                [[paste(sample_name,"_cell_meta",sep = "")]][["min_border_dist"]],
                                                row.names = colnames(tmp)))
  tmp <- AddMetaData(tmp, metadata = data.frame(dist_flag = ifelse(tmp@meta.data[["min_border_dist"]] < dist_cutoff * 3.3,"near","away"),
                                                row.names = colnames(tmp)))
  
  if(top_level!="all"){
    tmp_str <- as.character(top_level)
    tmp <- subset(tmp,subset = top_level == tmp_str)
  }
  tmp <- SetIdent(object = tmp, value  = tmp@meta.data[["dist_flag"]])
  marker_list <- FindMarkers(tmp,ident.1 = "near",ident.2 = "away",logfc.threshold = 0)
  marker_list <- marker_list[marker_list$p_val < p_cutoff,]
  #write.csv(marker_list,paste("near_plaque",sample_name,top_level,".csv",sep = "_"),row.names = T)
  return(marker_list)
}


ad_get_spatial_expr <- function(starmap_obj, 
                                sample_name, 
                                min_pct = 0.05, 
                                cell_meta, 
                                gene_list = NULL, logfc = NULL, logfc_cutoff = NULL,
                                filter = "logfc"){#or "he"
  min_pct <- 0.05
  starmap_obj <- starmap_raw
  cell_meta <- starmap@misc$AD_mouse9723_plaque_all$AD_mouse9723_cell_meta
  sample_name <- "AD_mouse9723"
  cell_meta <- starmap@misc$AD_mouse9494_plaque_all$AD_mouse9494_cell_meta
  sample_name <- "AD_mouse9494"
  
  starmap_obj <- starmap_64
  cell_meta <- starmap_64@misc[["AD_mouse9919_plaque_all"]][["AD_mouse9919_cell_meta"]]
  sample_name <- "AD_mouse9919"
  cell_meta <- starmap_64@misc[["AD_mouse9721_plaque_all"]][["AD_mouse9721_cell_meta"]]
  sample_name <- "AD_mouse9721"
  #Gene Filter if not provided
  #5% PCT
  if(is.null(gene_list)){
    gene_list <- row.names(starmap_obj@assays[["RNA"]]@counts)[
      rowSums(starmap_obj@assays[["RNA"]]@counts[,starmap_obj@meta.data[["sample"]] == sample_name] > 0)
      > sum(starmap_obj@meta.data[["sample"]] == sample_name) * min_pct]
  }

  tmp_df <- matrix(0, nrow = length(gene_list), ncol = 10) %>% as.data.frame(row.names = gene_list)
  colnames(tmp_df) <- paste("R",1:10,"0",sep = "")
  
  tmp_mat <- starmap_obj@assays[["RNA"]]@counts[gene_list,starmap_obj@meta.data[["sample"]] == sample_name]
  
  for (i in 1:10) {
    tmp_df[,i] <- rowSums(tmp_mat[,cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33]) /
      sum(cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33)
  }
  require(ComplexHeatmap)
  if(filter == "he"){
    he_gene_list <- gene_list[apply(tmp_df,1,FUN = function(x){
      vec <- as.numeric(x)
      return((vec[1]+vec[2])/2 > max(vec[3:10]))
    })]
  }else if(filter == "logfc"){
    he_gene_list <- gene_list[logfc > logfc_cutoff]
  }else if(filter == "sig_test"){
    tmp_df <- data.frame(tmp_df[,1:4], R50 = 0,sig_interval = 0, p_val = 1, gene = row.names(tmp_df),p_adj = 0)
    tmp_df$R50 <- rowSums(tmp_mat[,cell_meta$min_border_dist >= 4 * 33]) /
      sum(cell_meta$min_border_dist >= 4 * 33)
    tmp_df <- subset(tmp_df,apply(tmp_df[,1:4],1,max) > R50)
    tmp_df$sig_interval <- apply(tmp_df,MARGIN = 1,
                                 FUN = function(x){
                                   tmp_vec <- as.numeric(x[1:4])
                                   paste("R",which(tmp_vec == max(tmp_vec)),"0",sep = "")
                                 }) 
    tmp_df$p_val <- apply(tmp_df,MARGIN = 1,
                          FUN = function(x){
                            i <- switch(x[["sig_interval"]],
                                        "R10" = 1,"R20" = 1,"R30" = 3,"R40" = 4)
                            tmp_vec1 <- tmp_mat[x[["gene"]], cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33]
                            tmp_vec2 <- tmp_mat[x[["gene"]], !(cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33)]
                            return(wilcox.test(tmp_vec1,tmp_vec2)$p.value)
                          }) 
    for(i in 1:4){
      tmp_df$p_adj[tmp_df$sig_interval == paste("R",i,"0",sep = "")] <- 
        p.adjust(tmp_df$p_val[tmp_df$sig_interval == paste("R",i,"0",sep = "")])
    }
    tmp <- subset(tmp_df, p_val<0.01) %>% arrange(sig_interval,p_adj) %>%
      data.frame(pct = rowSums(tmp_mat[.$gene,])/dim(tmp_mat)[2])
    #For 64
    tmp <- subset(tmp_df, gene %in% tmp_vec_8) %>% arrange(sig_interval,p_adj)
    
    row.names(tmp) <- str_to_title(row.names(tmp))
    Heatmap(t(scale(t(tmp[,1:5]))),cluster_columns = F,cluster_rows  = F, show_row_names = T,column_title = paste(sample_name,min_pct)) 
    
    ####
    tmp_func <- function(cell_type, direction){
      gene_list <- switch(cell_type,
                          "Micro" = DE_13mo_Micro,
                          "Astro" = DE_13mo_Astro,
                          "Oligo" = DE_13mo_Oligo)
      if(direction == "UP"){
        gene_list <- gene_list %>% subset(subset = avg_logFC > 0) %>% .$rowname %>% .[1:10]
      }else{
        gene_list <- gene_list %>% subset(subset = avg_logFC < 0) %>% .$rowname %>% .[1:10]
      }
      tmp_df <- matrix(0, nrow = length(gene_list), ncol = 5) %>% as.data.frame(row.names = gene_list)
      colnames(tmp_df) <- paste("R",1:5,"0",sep = "")
      tmp_df <- data.frame(tmp_df,R60 = 0)
      
      tmp_mat <- starmap_obj@assays[["RNA"]]@counts[gene_list,starmap_obj@meta.data[["sample"]] == sample_name]
      
      for (i in 1:5) {
        tmp_df[,i] <- rowSums(tmp_mat[,cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33 & cell_meta$top_level == cell_type])  /
          sum(cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33 & cell_meta$top_level == cell_type)
      }
      tmp_df$R60 <- rowSums(tmp_mat[,cell_meta$min_border_dist >= 5 * 33 & cell_meta$top_level == cell_type]) /
        sum(cell_meta$min_border_dist >= 5 * 33 & cell_meta$top_level == cell_type)
      Heatmap(t(scale(t(tmp_df[,1:6]))),cluster_columns = F,cluster_rows  = T, show_row_names = T,column_title = paste(cell_type,direction))
    }
    
    
    
    pdf(file = "F3-5_tau_distr.pdf",width = 6,height = 4.5,useDingbats = FALSE)
    tmp_func("Micro","UP") %>% print()
    tmp_func("Astro","UP") %>% print()
    tmp_func("Oligo","UP") %>% print()
    tmp_func("Micro","DOWN") %>% print()
    tmp_func("Astro","DOWN") %>% print()
    tmp_func("Oligo","DOWN") %>% print()
    dev.off()
                                  
    
    ###                      
  }
  
  tmp_mat <- as.matrix(tmp_df[he_gene_list,1:10])
  overall_heatmap <- Heatmap(log2(1+as.matrix(tmp_df[,1:10])),cluster_columns = F, show_row_names = F ) 
  raw_he_heatmap <- Heatmap(tmp_mat[,1:5],cluster_columns = F, row_title = sample_name)
  scaled_he_heatmap <- Heatmap(t(scale(t(tmp_mat[,1:5]))),cluster_columns = F, row_title = sample_name)
  return(list("spatial_expr" = tmp_df,
              "gene_list" = gene_list,
              "he_gene_list" = he_gene_list,
              "overall_heatmap" = overall_heatmap,
              "raw_he_heatmap" = raw_he_heatmap,
              "scaled_he_heatmap" = scaled_he_heatmap))
  
}

#####Figure Export#####
#ad_distr_fig_export(starmap,"AD_mouse9494",cell_type_list,cell_type_palette_list)
#ad_distr_fig_export(starmap,"AD_mouse9723",cell_type_list,cell_type_palette_list,scale_max = c("Micro"=1,"Astro"=0.25,"Oligo"=0.75, "Ex"=0.6,"Inhi"=0.15))
ad_distr_fig_export_all <- function(starmap_obj, 
                                    sample_name, 
                                    cell_type_list, 
                                    cell_type_palette_list, 
                                    scale_max = c("Micro"=1,"Astro"=0.2,"Oligo"=0.4, "Ex"=0.6,"Inhi"=0.15)){#M,A,O,Ex,Inhi
  pdf(file = paste(sample_name,"cell_distr_fig.pdf",sep = "_"),width = 12,height = 9)
  plot_list <- list()
  #Top Level
  print(starmap_obj@misc[[paste(sample_name,"plaque_all",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[3]] + 
          labs(title = paste(sample_name,"all Percentage")))
  
  #Sub Cluster
  for(i in c("Micro","Astro","Ex","Inhi")){
    palette_list <- cell_type_palette_list[[i]]
    names(palette_list) <- cell_type_list[[i]]
    print(ad_hl_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_all",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                           starmap_obj@misc[[paste(sample_name,"plaque_all",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                           palette_list = palette_list,
                           scale = "percentage",title = paste(sample_name,i,"all"),y_min = 0,y_max = scale_max[[i]]))
  }
  
  #OPC/Oligo
  palette_list <- cell_type_palette_list$Oligo
  names(palette_list) <- cell_type_list$Oligo
  tmp_df <- data.frame(min_border_dist = starmap_obj@misc[[paste(sample_name,"plaque_all",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                       cell_type =       starmap_obj@misc[[paste(sample_name,"plaque_all",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type)
  tmp_df$cell_type[tmp_df$cell_type %in% cell_type_list$OPC] <- "OPC" 
  print(ad_hl_cell_distr(tmp_df$min_border_dist,
                         tmp_df$cell_type,
                         palette_list = palette_list,
                         scale = "percentage",title = paste(sample_name,"Oligo/OPC","all"),y_min = 0,y_max = scale_max[["Oligo"]]))
  dev.off()
}

ad_distr_fig_export_regional <- function(starmap_obj, 
                                         sample_name, 
                                         cell_type_list, 
                                         cell_type_palette_list, 
                                         scale_max = c("Micro"=1,"Astro"=0.2,"Oligo"=0.9, "Ex"=0.6,"Inhi"=0.15, "CADG" = 0.6),
                                         scale_density_max = NULL){#M,A,O,Ex,Inhi 
  pdf(file = paste(sample_name,"cell_distr_fig.pdf",sep = "_"),width = 12,height = 9)
  plot_list <- list()
  #Top Level
  print(starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[3]] + 
          labs(title = paste(sample_name,"Cortex Percentage")))
  print(starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[3]] + 
          labs(title = paste(sample_name,"Hippocampus Percentage")))
  print(starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[3]] + 
          labs(title = paste(sample_name,"White Matter Percentage")))
  
  print(starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[2]] + 
          labs(title = paste(sample_name,"Cortex Density")))
  print(starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[2]] + 
          labs(title = paste(sample_name,"Hippocampus Density")))
  print(starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"t_l_d_plotlist",sep="_")]][[2]] + 
          labs(title = paste(sample_name,"White Matter Density")))
  dev.off()
  tmp_get_img_size <- function(region, region_img){
    img_size <- 0
    if(region == "all"){
      img_size <- dim(region_img)[1] * dim(region_img)[2]
    }else if(region == "Sub-cortical"){
      img_size <- sum(region_img == 2) + 
        sum(region_img == 3)
    }else{
      x <- switch(region,
                  "Cortex" = 1,"White Matter" = 2, "Hippocampus" = 3)
      img_size <- sum(region_img == x)
    }
    
    return(img_size)
  }
  
  
  region_img <- starmap_obj@misc[[paste(sample_name,"morph",sep = "_")]][["region"]]
  area_vec_c = starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"plaque_dilate_area",sep = "_")]]
  area_vec_wm = starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"plaque_dilate_area",sep = "_")]]
  area_vec_h = starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"plaque_dilate_area",sep = "_")]]
  #Sub Cluster
  for(i in c("Micro","Astro","Ex","Inhi","CADG")){
    pdf(file = paste(sample_name,i,"cell_distr_fig.pdf",sep = "_"),width = 12,height = 9)
    palette_list <- cell_type_palette_list[[i]]
    names(palette_list) <- cell_type_list[[i]]
    print(ad_hl_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                           starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                           palette_list = palette_list,
                           scale = "percentage",title = paste(sample_name,i,"Cortex Percentage"),y_min = 0,y_max = scale_max[[i]]))
    print(ad_plot_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                             starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                             palette_list = palette_list, area_vec = area_vec_c,img_size = tmp_get_img_size("Cortex",region_img),
                             scale = "density",title = paste(sample_name,i,"Cortex Density")))
    
    print(ad_hl_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                           starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                           palette_list = palette_list,
                           scale = "percentage",title = paste(sample_name,i,"White Matter Percentage"),y_min = 0,y_max = scale_max[[i]]))
    print(ad_plot_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                             starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                             palette_list = palette_list, area_vec = area_vec_wm,img_size = tmp_get_img_size("White Matter",region_img),
                             scale = "density",title = paste(sample_name,i,"White Matter Density")))
    
    print(ad_hl_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                           starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                           palette_list = palette_list,
                           scale = "percentage",title = paste(sample_name,i,"Hippocampus Percentage"),y_min = 0,y_max = scale_max[[i]]))
    print(ad_plot_cell_distr(starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                             starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type,
                             palette_list = palette_list, area_vec = area_vec_h,img_size = tmp_get_img_size("Hippocampus",region_img),
                             scale = "density",title = paste(sample_name,i,"Hippocampus Density")))
    
    dev.off()
  }
  
  #OPC/Oligo
  pdf(file = paste(sample_name,"OPC_n_Oligo","cell_distr_fig.pdf",sep = "_"),width = 12,height = 9)
  palette_list <- cell_type_palette_list$Oligo
  names(palette_list) <- cell_type_list$Oligo
  tmp_df <- data.frame(min_border_dist = starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                       cell_type =       starmap_obj@misc[[paste(sample_name,"plaque_Cortex",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type)
  tmp_df$cell_type[tmp_df$cell_type %in% cell_type_list$OPC] <- "OPC" 
  print(ad_hl_cell_distr(tmp_df$min_border_dist,
                         tmp_df$cell_type,
                         palette_list = palette_list,
                         scale = "percentage",title = paste(sample_name,"Oligo/OPC","Cortex Percentage"),y_min = 0,y_max = scale_max[["Oligo"]]))
  print(ad_plot_cell_distr(tmp_df$min_border_dist,
                           tmp_df$cell_type,
                           palette_list = palette_list, area_vec = area_vec_c,img_size = tmp_get_img_size("Cortex",region_img),
                           scale = "density",title = paste(sample_name,"Oligo/OPC","Cortex Density")))
  
  tmp_df <- data.frame(min_border_dist = starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                       cell_type =       starmap_obj@misc[[paste(sample_name,"plaque_White Matter",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type)
  tmp_df$cell_type[tmp_df$cell_type %in% cell_type_list$OPC] <- "OPC" 
  print(ad_hl_cell_distr(tmp_df$min_border_dist,
                         tmp_df$cell_type,
                         palette_list = palette_list,
                         scale = "percentage",title = paste(sample_name,"Oligo/OPC","White Matter Percentage"), y_min = 0,y_max = scale_max[["Oligo"]]))
  print(ad_plot_cell_distr(tmp_df$min_border_dist,
                           tmp_df$cell_type,
                           palette_list = palette_list, area_vec = area_vec_wm,img_size = tmp_get_img_size("White Matter",region_img),
                           scale = "density",title = paste(sample_name,"Oligo/OPC","White Matter Density")))
  
  tmp_df <- data.frame(min_border_dist = starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$min_border_dist,
                       cell_type =       starmap_obj@misc[[paste(sample_name,"plaque_Hippocampus",sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]]$cell_type)
  tmp_df$cell_type[tmp_df$cell_type %in% cell_type_list$OPC] <- "OPC" 
  print(ad_hl_cell_distr(tmp_df$min_border_dist,
                         tmp_df$cell_type,
                         palette_list = palette_list,
                         scale = "percentage",title = paste(sample_name,"Oligo/OPC","Hippocampus Percentage"), y_min = 0,y_max = scale_max[["Oligo"]]))
  print(ad_plot_cell_distr(tmp_df$min_border_dist,
                           tmp_df$cell_type,
                           palette_list = palette_list, area_vec = area_vec_h,img_size = tmp_get_img_size("Hippocampus",region_img),
                           scale = "density",title = paste(sample_name,"Oligo/OPC","Hippocampus Density")))
  dev.off()
}


ad_calc_cluster_dist <- function(cell_meta,
                                 sample_name,
                                 region = "all",
                                 Top_level_clusters = NULL,
                                 Sub_level_clusters = NULL,
                                 plaque_meta = NULL,
                                 shuffle = F){#Top_level_cluster should be a vector contain cluster names
  cell_meta <- cell_meta[cell_meta[["sample"]] == sample_name,]
  if(region == "Cortex"){
    cell_meta <- cell_meta[cell_meta$region == "Cortex",]
    plaque_meta <- plaque_meta[plaque_meta$region == "Cortex",]
  }else if(region == "White Matter"){
    cell_meta <- cell_meta[cell_meta$region == "White Matter",]
    plaque_meta <- plaque_meta[plaque_meta$region == "White Matter",]
  }else if(region == "Hippocampus"){
    cell_meta <- cell_meta[cell_meta$region == "Hippocampus",]
    plaque_meta <- plaque_meta[plaque_meta$region == "Hippocampus",]
  }else if(region == "Sub-cortical"){
    cell_meta <- cell_meta[cell_meta$region %in% c("Hippocampus","White Matter"),]
    plaque_meta <- plaque_meta[plaque_meta$region %in% c("Hippocampus","White Matter"),]
  }
  if(!is.null(Top_level_clusters)){
    if(is.null(plaque_meta)){
      cell_meta <- data.frame(x = c(cell_meta$x),
                              y = c(cell_meta$y),
                              top_level = c(as.character(cell_meta$top_level)))
      Top_level_clusters <- c(as.character(Top_level_clusters))
    }else{
      cell_meta <- data.frame(x = c(cell_meta$x,plaque_meta$m.cx),
                              y = c(cell_meta$y,plaque_meta$m.cy),
                              top_level = c(as.character(cell_meta$top_level),rep("plaque",dim(plaque_meta)[1])))
      Top_level_clusters <- c(as.character(Top_level_clusters),"plaque")
    }
    if(shuffle){
      set.seed(123)
      cell_meta$top_level[cell_meta$top_level!="plaque"] <- sample(cell_meta$top_level[cell_meta$top_level!="plaque"])
    }
    
    for(i in 1:length(Top_level_clusters)){
      cell_meta <- cbind(cell_meta, tmp = -1)
      colnames(cell_meta)[dim(cell_meta)[2]] = as.character(Top_level_clusters[i])
    }
    
    #Calc Min Dist
    print("Calc Min cluster Dist")
    for(i in 1:length(Top_level_clusters)){
      print(paste("Calculating:", Top_level_clusters[i],", Progress:",i,'/',length(Top_level_clusters)))
      cell_meta[[Top_level_clusters[i]]] <- apply(cell_meta,
                                                  MARGIN = 1,
                                                  FUN = function(x){
                                                    if(as.character(x[3]) == as.character(Top_level_clusters[i])){
                                                      tmp <- sqrt((as.numeric(x[1]) - 
                                                                     cell_meta$x[as.character(cell_meta$top_level) == 
                                                                                   as.character(Top_level_clusters[i])])^2+
                                                                    (as.numeric(x[2]) - 
                                                                       cell_meta$y[as.character(cell_meta$top_level) == 
                                                                                     as.character(Top_level_clusters[i])])^2)
                                                      return(min(tmp[tmp > 0.0001]))
                                                    }else{
                                                      tmp <- sqrt((as.numeric(x[1]) - 
                                                                     cell_meta$x[as.character(cell_meta$top_level) == 
                                                                                   as.character(Top_level_clusters[i])])^2+
                                                                    (as.numeric(x[2]) - 
                                                                       cell_meta$y[as.character(cell_meta$top_level) == 
                                                                                     as.character(Top_level_clusters[i])])^2)
                                                      return(min(tmp))
                                                    }
                                                  })
    }
    dist_stat <- data.frame(type_from = rep(Top_level_clusters,each = length(Top_level_clusters)),
                            type_to = rep(Top_level_clusters,times = length(Top_level_clusters)),
                            mean_dist = 0)
    dist_stat$mean_dist <- apply(dist_stat,1,
                                 FUN = function(x){
                                   return(mean(cell_meta[as.character(cell_meta$top_level) == as.character(x[1]),as.character(x[2])]))
                                 })
  }else{
    if(is.null(plaque_meta)){
      cell_meta <- data.frame(x = c(cell_meta$x),
                              y = c(cell_meta$y),
                              cell_type = c(as.character(cell_meta$cell_type)))
      Sub_level_clusters <- c(as.character(Sub_level_clusters))
    }else{
      cell_meta <- data.frame(x = c(cell_meta$x,plaque_meta$m.cx),
                              y = c(cell_meta$y,plaque_meta$m.cy),
                              cell_type = c(as.character(cell_meta$cell_type),rep("plaque",dim(plaque_meta)[1])))
      Sub_level_clusters <- c(as.character(Sub_level_clusters),"plaque")
    }
    Sub_level_clusters <- as.character(Sub_level_clusters)
    for(i in 1:length(Sub_level_clusters)){
      cell_meta <- cbind(cell_meta, tmp = -1)
      colnames(cell_meta)[dim(cell_meta)[2]] = as.character(Sub_level_clusters[i])
    }
    if(shuffle){
      set.seed(123)
      cell_meta$cell_type <- sample(cell_meta$cell_type)
    }
    #Calc Min Dist
    print("Calc Min cluster Dist")
    for(i in 1:length(Sub_level_clusters)){
      print(paste("Calculating:", Sub_level_clusters[i],", Progress:",i,'/',length(Sub_level_clusters)))
      cell_meta[[Sub_level_clusters[i]]] <- apply(cell_meta,
                                                  MARGIN = 1,
                                                  FUN = function(x){
                                                    if(as.character(x[["cell_type"]]) == as.character(Sub_level_clusters[i])){
                                                      tmp <- sqrt((as.numeric(x[["x"]]) - 
                                                                     cell_meta$x[as.character(cell_meta$cell_type) == 
                                                                                   as.character(Sub_level_clusters[i])])^2+
                                                                    (as.numeric(x[["y"]]) - 
                                                                       cell_meta$y[as.character(cell_meta$cell_type) == 
                                                                                     as.character(Sub_level_clusters[i])])^2)
                                                      return(min(tmp[tmp > 0.0001]))
                                                    }else{
                                                      tmp <- sqrt((as.numeric(x[["x"]]) - 
                                                                     cell_meta$x[as.character(cell_meta$cell_type) == 
                                                                                   as.character(Sub_level_clusters[i])])^2+
                                                                    (as.numeric(x[["y"]]) - 
                                                                       cell_meta$y[as.character(cell_meta$cell_type) == 
                                                                                     as.character(Sub_level_clusters[i])])^2)
                                                      return(min(tmp))
                                                    }
                                                  })
    }
    dist_stat <- data.frame(type_from = rep(Sub_level_clusters,each = length(Sub_level_clusters)),
                            type_to = rep(Sub_level_clusters,times = length(Sub_level_clusters)),
                            mean_dist = 0)
    dist_stat$mean_dist <- apply(dist_stat,1,
                                 FUN = function(x){
                                   return(mean(cell_meta[as.character(cell_meta$cell_type) == as.character(x[1]),as.character(x[2])]))
                                 })
    if(shuffle){
      for(i in 1:length(Sub_level_clusters)){
        dist_stat$mean_dist[dist_stat$type_to == Sub_level_clusters[i]] <- 
          mean(dist_stat$mean_dist[dist_stat$type_to == Sub_level_clusters[i]])
      }
    }
  }
  dist_stat$mean_dist <- dist_stat$mean_dist/3.3
  #3.3pixel = 1um
  
  return(dist_stat)
}
ad_plot_type_dist <- function(disease_dist, ctrl_dist,disease_name,ctrl_name){
  require(RColorBrewer)
  plot_list <- list()
  dist_limit <- max(c(disease_dist$mean_dist,ctrl_dist$mean_dist))
  plot_list[[1]] <-  ggplot(data = disease_dist,aes(x = type_from, y=type_to, fill = mean_dist)) +
    geom_tile()+
    labs(title = disease_name) +
    scale_fill_gradientn(colours = brewer.pal(10,"RdYlGn"), 
                         limit = c(0,dist_limit), space = "Lab", #max(log2(1+AD9723_cluster_dist_stat$mean_dist))
                         name="Distance(micron)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 8, hjust = 1),
          legend.position = "top")
  plot_list[[2]] <-  ggplot(data = ctrl_dist,aes(x = type_from, y=type_to, fill = mean_dist)) +
    geom_tile()+
    labs(title = ctrl_name) +
    scale_fill_gradientn(colours = brewer.pal(10,"RdYlGn"),
                         limit = c(0,dist_limit), space = "Lab", 
                         name="Distance(micron)") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 8, hjust = 1),
          legend.position = "top")
  
  return(ggarrange(plotlist = plot_list,nrow = 1,ncol = 2))
  
}


