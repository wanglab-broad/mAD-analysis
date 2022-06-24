#####Function Definition: preprocess and calc cell interval#####
#TEST
#starmap_obj       <- starmap_raw
#img_list          <- starmap_obj@misc$ADmouse_9494_morph
#starmap_obj@misc  <- list()
#sample_name       <- "ADmouse_9494"
#min_area          <- 400
#region            <- "Cortex"
#palette_list <- as.character(unlist(starmap_raw@misc[["top_hex_dict"]]))[colnames(starmap_raw@misc[["top_hex_dict"]]) != "DG"]
#names(palette_list) <- colnames(starmap_raw@misc[["top_hex_dict"]])[colnames(starmap_raw@misc[["top_hex_dict"]]) != "DG"]
#rm(list = c("starmap_obj","img_list","cell_meta","plaque_meta","sample_name","min_area","region"))        
#
ad_calc_plaque_dist2<- function(starmap_obj,
                                img_list,
                                min_area = 400,
                                sample_name,
                                palette_list,
                                region = "all"){
  #Calc Plaque Center
  require(EBImage)
  message("Load Image")
  plaque <- Image(img_list[["plaque"]],colormode = Grayscale) %>% bwlabel()
  message("Compute Plaque Features")
  plaque_meta <- cbind(computeFeatures.moment(plaque),computeFeatures.shape(plaque)) %>% 
    as.data.frame()  %>% data.frame(., plaque_id = as.numeric(row.names(.)))
  
  cell_meta <- data.frame(starmap_obj@meta.data[starmap_obj@meta.data[["sample"]] == sample_name,],nearest_plaque = 0)
  
  #Organize cell_meta
  #orig_index |x  |y  |batch  |time |group  |top_level  |cell_type  |region |min_dist |
  if(is.null(cell_meta$region)){
    cell_meta <-data.frame(orig_index = cell_meta$orig_index,
                           x = cell_meta$x, y = cell_meta$y,
                           batch = cell_meta$batch, time = cell_meta$time, group = cell_meta$group,
                           top_level = cell_meta$top_level, cell_type = cell_meta$cell_type,
                           region = NA,
                           nearest_plaque = 0, interval = "other")
  }else{
    cell_meta <-data.frame(orig_index = cell_meta$orig_index,
                           x = cell_meta$x, y = cell_meta$y,
                           batch = cell_meta$batch, time = cell_meta$time, group = cell_meta$group,
                           top_level = cell_meta$top_level, cell_type = cell_meta$cell_type,
                           region = cell_meta$region,
                           nearest_plaque = 0, interval = "other")
  }
  
  #Region Selection
  if(region != "all"){
    plaque_meta <- data.frame(plaque_meta,region = -1)
    plaque_meta$region <- apply(plaque_meta,1,
                                FUN = function(x){
                                  switch(img_list[["region"]][round(as.numeric(x[1])),round(as.numeric(x[2]))],
                                         '1' = "Cortex", '2' = "Corpus Callosum", '3' = "Hippocampus")
                                })
    if(sum(is.na(cell_meta$region)==F) == 0){
      cell_meta$region <- apply(cell_meta,1,
                                FUN = function(x){
                                  switch(img_list[["region"]][round(as.numeric(x[2])),round(as.numeric(x[3]))],
                                         '1' = "Cortex", '2' = "Corpus Callosum", '3' = "Hippocampus")
                                })
    }
  }
  #Save Raw palque meta for dilate
  plaque_meta_bk <- plaque_meta
  plaque_meta <- plaque_meta[plaque_meta$s.area > min_area,]
  
  if(region == "Sub-cortical"){
    cell_meta <- cell_meta[cell_meta$region %in% c("Hippocampus","Corpus Callosum"),]
    plaque_meta <- plaque_meta[plaque_meta$region %in% c("Hippocampus","Corpus Callosum"),]
  }else if(region != "all"){
    cell_meta <- cell_meta[cell_meta$region == region,]
    plaque_meta <- plaque_meta[plaque_meta$region == region,]
  }
  
  
  #Calc Plaque 5-round dilate area 
  message("Calc size of dilated areas")
  if(region == "Sub-cortical"){#Filter out 
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area | 
                                                          !(plaque_meta_bk$region %in% c("Hippocampus","Corpus Callosum"))],
                        reenumerate = F)
  }else if(region != "all"){
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area | 
                                                          plaque_meta_bk$region != region], reenumerate = F)
  }else{
    plaque <- rmObjects(plaque,plaque_meta_bk$plaque_id[plaque_meta_bk$s.area <= min_area], reenumerate = F)
  }
  
  #Do voronoi
  voronoi_img <- propagate(seeds = plaque, x = Image(dim = dim(plaque)))
  
  cell_meta$nearest_plaque <- apply(cell_meta,1,
                                    FUN = function(x){
                                      as.numeric(voronoi_img[round(as.numeric(x[["x"]])),round(as.numeric(x[["y"]]))])
                                    })
  
  if(region == "all"){
    img_size <- dim(plaque)[1] * dim(plaque)[2]
  }else if(region == "Sub-cortical"){
    img_size <- sum(img_list[["region"]] == 2) + 
      sum(img_list[["region"]] == 3)
    region_img <- Image(img_list[["region"]] == 2 |img_list[["region"]] == 3)
  }else{
    x <- switch(region,
                "Cortex" = 1,"Corpus Callosum" = 2, "Hippocampus" = 3)
    img_size <- sum(img_list[["region"]] == x)
    region_img <- Image(img_list[["region"]] == x)
  }
  
  #Assign interval
  dilate_ls <- list()
  size_vec  <- rep(0,5); names(size_vec) <- paste(1:5,"0um",sep = "")
  dilate_img <- plaque > 0 
  dilate_img_bk <- dilate_img
  for(i in 1:5){
    dilate_img_bk <- dilate_img
    dilate_img <- dilate(plaque > 0 , kern = makeBrush(round(63.5*i), shape = "disc"))
    if(region != "all") dilate_img <- dilate_img * region_img
    if(i != 1) dilate_img2 <- dilate_img * (1 - dilate_img_bk) else dilate_img2 <- dilate_img
    size_vec[i] <- sum(dilate_img2) 
    
    cell_meta[["interval"]] <- apply(cell_meta,1,
                                     FUN = function(x){
                                       if(dilate_img2[as.numeric(x[["x"]]),as.numeric(x[["y"]])] != 0 ) return(paste(i,"0um",sep = ""))
                                       else return(x[["interval"]])
                                     })
    dilate_img2 <- dilate_img2 * voronoi_img
    dilate_ls[[paste0(i,"0um_dilate_info")]] <- as.data.frame(computeFeatures.shape(dilate_img2)) %>% 
      data.frame(plaque_id = as.numeric(row.names(.)), .)
  }
  
  #size_vec[1] <- tmp_vec[1] - sum(plaque)
  
  return_list <- list()
  return_list[[paste(sample_name,"plaque_dilate_area",sep = "_")]] <- size_vec
  return_list[[paste(sample_name,"cell_meta",sep = "_")]] <- cell_meta
  return_list[[paste(sample_name,"plaque_meta",sep = "_")]] <- plaque_meta
  return_list[[paste(sample_name,"plaque_img",sep = "_")]] <- plaque
  return_list[[paste(sample_name,"plaque_dilate_info",sep = "_")]] <- dilate_ls


  starmap_obj@misc[[paste(sample_name,"plaque",region,sep = "_")]] <- return_list
  return(starmap_obj)
}
