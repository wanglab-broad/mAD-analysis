#####Function Definition: Statistical Analysis#####
######Test Para
#res_ls <- starmap_raw@misc$ADmouse_9494_plaque_Cortex
#sample_name <- "ADmouse_9494"
#region_size <- sum(starmap_raw@misc$ADmouse_9494_morph$region == 3)
#region <- "Hippocampus"
#top_level_list <- c("CTX-Ex", "Endo", "Inh", "Micro", "Astro", "OPC","Oligo", "SMC") 
#
#res_ls <- starmap_64@misc$AD_mouse9721_plaque_all
#sample_name <- "AD_mouse9721"
#region_size <-  dim(starmap_64@misc[[paste0(sample_name,"_morph")]]$plaque)[1] *
#  dim(starmap_64@misc[[paste0(sample_name,"_morph")]]$plaque)[2]
#region <- "all"
#top_level_list <- c("CTX-Ex", "CA1","CA2","CA3", "DG","Endo", "Inh", "Micro", "Astro", "OPC","Oligo", "SMC")
######Func1
ad_distr_stat <- function(res_ls, sample_name, region_size, cell_type_list, top_level_list, region = "all"){
  dilate_ls <- res_ls[[paste0(sample_name,"_plaque_dilate_info")]]
  cell_meta <- res_ls[[paste0(sample_name,"_cell_meta")]]
  plaque_meta <- res_ls[[paste0(sample_name,"_plaque_meta")]]
  dilate_area <- res_ls[[paste0(sample_name,"_plaque_dilate_area")]]
  cell_meta <- subset(cell_meta, top_level %in% top_level_list)
  return_ls <- list()
  return_ls[["avg_stat"]] <- list()
  #Top level
  tmp_vec <- setNames(sapply(top_level_list,function(x){1000000* sum(cell_meta$top_level == x)/(region_size/res_factor^2)}), 
                      top_level_list)#avg density
  
  tmp_df <- data.frame(plaque_id = rep(plaque_meta$plaque_id, each = length(top_level_list)),
                       top_level = rep(top_level_list, nrow(plaque_meta)))
  for(i in 1:5){
    tmp_df[[paste0("R",i*10)]] <- apply(tmp_df,1,
                                        FUN = function(x){
                                          if(!as.numeric(x[["plaque_id"]]) %in%
                                             dilate_ls[[paste0(i,"0um_dilate_info")]]$plaque_id){ NA
                                          }else{
                                            1000000 * sum(cell_meta$interval == paste0(i,"0um") & 
                                                            cell_meta$nearest_plaque == as.numeric(x[["plaque_id"]]) &
                                                            cell_meta$top_level == as.character(x[["top_level"]])) /                                            (dilate_ls[[paste0(i,"0um_dilate_info")]]$s.area[dilate_ls[[paste0(i,"0um_dilate_info")]]$plaque_id ==  as.numeric(x[["plaque_id"]])]/res_factor^2)
                                          }
                                        })
    
  }
  return_ls[["top_level"]] <- tmp_df
  return_ls[["avg_stat"]][["top_level"]] <- tmp_vec
  #Sub type
  for(i in names(cell_type_list)){
    tmp_df <- data.frame(plaque_id = rep(plaque_meta$plaque_id, each = length(cell_type_list[[i]])),
                         subtype = rep(cell_type_list[[i]], nrow(plaque_meta)))
    for(j in 1:5){
      tmp_df[[paste0("R",j*10)]] <- apply(tmp_df,1,
                                          FUN = function(x){
                                            if(!as.numeric(x[["plaque_id"]]) %in%
                                               dilate_ls[[paste0(j,"0um_dilate_info")]]$plaque_id) NA
                                            else{
                                              1000000 * sum(cell_meta$interval == paste0(j,"0um") & 
                                                              cell_meta$nearest_plaque == as.numeric(x[["plaque_id"]]) &
                                                              cell_meta$cell_type == as.character(x[["subtype"]])) /                                          (dilate_ls[[paste0(j,"0um_dilate_info")]]$s.area[dilate_ls[[paste0(j,"0um_dilate_info")]]$plaque_id ==  as.numeric(x[["plaque_id"]])]/res_factor^2)
                                            }
                                          }) %>% unlist
    }
    
    tmp_vec <- setNames(sapply(cell_type_list[[i]],function(x){1000000* sum(cell_meta$cell_type == x)/(region_size/res_factor^2)}), 
                        cell_type_list[[i]])#avg density
    return_ls[[i]] <- tmp_df
    return_ls[["avg_stat"]][[i]] <- tmp_vec
  }
  return_ls[["cell_meta"]] <- cell_meta
  return_ls[["dilate_vec"]] <- dilate_area
  return_ls
}
######Test Para
#sample_list <- c("ADmouse_9494","ADmouse_11346") 
#res_ls <- distr_stat_ls
#region = "Hippocampus"
#
#sample_list <- c("AD_mouse9721") 
#res_ls <- distr_stat_ls
#region = "all"
#####Func2
ad_distr_combined_analysis <- function(sample_list, res_ls, region, cell_type_list){
  return_ls <- list()
  return_ls[["plot_ls"]] <- list()
  return_ls[["stat_ls"]] <- list()
  stacked_vln <- function(df, palette_list, title, level, mean_vec){
    colnames(df)[colnames(df) == level] <- "cell_type"
    ggplot(df,aes(variable, value, fill = cell_type)) +
      geom_violin(scale = "width", adjust = 1, trim = TRUE) +
      scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
        c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
      facet_grid(rows = vars(cell_type) , scales = "free", switch = "y") +
      geom_hline(data = data.frame(cell_type = factor(names(mean_vec)), value = tmp_vec), 
                 aes(yintercept = value),color = "green") + 
      theme_cowplot(font_size = 12) +
      theme(legend.position = "none", panel.spacing = unit(0, "lines"),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = palette_list,
                        breaks = names(palette_list)) + 
      ggtitle(paste(paste(sample_list,collapse = " "), title)) + xlab("Interval") + ylab("Density")
  }
  
  cell_meta <- res_ls[[paste(sample_list[1],region,sep = "_")]][["cell_meta"]]
  dilate_vec <- res_ls[[paste(sample_list[1],region,sep = "_")]][["dilate_vec"]]
  tmp_df <- as.data.frame(res_ls[[paste(sample_list[1],region,sep = "_")]][["top_level"]])
  tmp_vec <- res_ls[[paste(sample_list[1],region,sep = "_")]][["avg_stat"]][["top_level"]]
  if(length(sample_list) > 1){
    for(i in 2:length(sample_list)){
      cell_meta <- rbind(cell_meta,res_ls[[paste(sample_list[i],region,sep = "_")]][["cell_meta"]])
    }
    
    for(i in 2:length(sample_list)){
      dilate_vec <- dilate_vec + res_ls[[paste(sample_list[i],region,sep = "_")]][["dilate_vec"]]
    }
    #Top level
    
    for(i in 2:length(sample_list)){
      tmp_df <- rbind(tmp_df,as.data.frame(res_ls[[paste(sample_list[i],region,sep = "_")]][["top_level"]]))
    }
    
    for(i in 2:length(sample_list)){
      tmp_vec <- tmp_vec+res_ls[[paste(sample_list[i],region,sep = "_")]][["avg_stat"]][["top_level"]]
    }
  }
  
  tmp_vec <- tmp_vec/length(sample_list)
  tmp_list <- list()
  tmp_list[["perc_df"]] <- data.frame(top_level = levels(as.factor(tmp_df$top_level)))
  tmp_list[["dens_df"]] <- data.frame(top_level = levels(as.factor(tmp_df$top_level)))
  tmp_list[["stat_c_df"]] <- data.frame(top_level = levels(as.factor(tmp_df$top_level)))#Chi-sq test
  tmp_list[["stat_t_df"]] <- data.frame(top_level = levels(as.factor(tmp_df$top_level)))#one sample t-test
  for(i in 1:5){
    tmp_list[["perc_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["perc_df"]],1,
                                                       FUN = function(x){
                                                         sum(cell_meta$top_level == x[["top_level"]] &
                                                               cell_meta$interval == paste0(i,"0um"))/
                                                           sum(cell_meta$interval == paste0(i,"0um"))})
    tmp_list[["dens_df"]][[paste0("R",i*10)]]  <- apply(tmp_list[["dens_df"]],1,
                                                        FUN = function(x){
                                                          1000000 * sum(cell_meta$top_level == x[["top_level"]] &
                                                                          cell_meta$interval == paste0(i,"0um"))/
                                                            (dilate_vec[[paste0(i,"0um")]]/res_factor^2) })
    
    tmp_list[["stat_t_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["stat_t_df"]],1,
                                                         FUN = function(x){
                                                           tmp_vec2 <- as.numeric(tmp_df[[paste0("R",i*10)]]) %>% 
                                                             .[tmp_df$top_level == x[["top_level"]]] %>% na.omit()
                                                           
                                                           tmp_obj <- t.test(x = tmp_vec2,
                                                                             mu = tmp_vec[[x[["top_level"]]]], 
                                                                             alternative = ifelse(mean(tmp_vec2) > tmp_vec[[x[["top_level"]]]],
                                                                                                  "greater", "less"))
                                                           tmp_obj$p.value
                                                         })
    
    tmp_list[["stat_c_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["stat_c_df"]],1,
                                                         FUN = function(x){
                                                           a = sum(cell_meta$top_level == x[["top_level"]] & cell_meta$interval == paste0(i,"0um"))
                                                           b = sum(cell_meta$top_level != x[["top_level"]] & cell_meta$interval == paste0(i,"0um"))
                                                           c = sum(cell_meta$top_level == x[["top_level"]] & cell_meta$interval != paste0(i,"0um"))
                                                           d = sum(cell_meta$top_level != x[["top_level"]] & cell_meta$interval != paste0(i,"0um"))
                                                           
                                                           tmp_obj <- suppressWarnings(chisq.test(x = matrix(c(a,c,b,d), nrow = 2)))
                                                           tmp_obj$p.value
                                                         })
  }
  #Add overall 
  tmp_list[["dens_df"]] <- left_join(tmp_list[["dens_df"]],
                                     data.frame(top_level = names(tmp_vec),overall = tmp_vec), by = "top_level")
  tmp_list[["perc_df"]][["overall"]] <- apply(tmp_list[["perc_df"]],1,
                                              FUN = function(x){
                                                sum(cell_meta$top_level == x[["top_level"]])/nrow(cell_meta)})
  return_ls[["stat_ls"]][["combined"]] <- tmp_list
  #Stacked Barplot
  require(reshape2)
  require(cowplot)
  require(magrittr)
  tmp_df %<>% melt(id.vars = c("plaque_id","top_level")) %>% subset(!is.na(value))
  for(i in levels(factor(tmp_df$top_level))){
    if(sum(tmp_df$value == 0 & tmp_df$top_level == i) / sum(tmp_df$value == 0 & tmp_df$top_level == i) < 0.99){
      tmp_df <- tmp_df[tmp_df$top_level != i|tmp_df$value < quantile(tmp_df$value[tmp_df$top_level == i], probs = 0.99),]
    }
  }
  return_ls[["plot_ls"]][["combined"]] <- stacked_vln(tmp_df, top_level_palette_list, "Top Level", "top_level", tmp_vec)
  
  for(j in names(cell_type_list)){
    
    tmp_df <- as.data.frame(res_ls[[paste(sample_list[1],region,sep = "_")]][[j]])
    tmp_vec <- res_ls[[paste(sample_list[1],region,sep = "_")]][["avg_stat"]][[j]]
    if(length(sample_list) > 1){
      for(i in 2:length(sample_list)){
        tmp_df <- rbind(tmp_df,as.data.frame(res_ls[[paste(sample_list[i],region,sep = "_")]][[j]]))
      }
      
      for(i in 2:length(sample_list)){
        tmp_vec <- tmp_vec+res_ls[[paste(sample_list[i],region,sep = "_")]][["avg_stat"]][[j]]
      }
    }
    
    tmp_vec <- tmp_vec/length(sample_list)
    tmp_list <- list()
    tmp_list[["perc_df"]]  <- data.frame(subtype = levels(as.factor(tmp_df$subtype)))
    tmp_list[["dens_df"]]  <- data.frame(subtype = levels(as.factor(tmp_df$subtype)))
    tmp_list[["stat_c_df"]] <- data.frame(subtype = levels(as.factor(tmp_df$subtype)))#Chi-sq test
    tmp_list[["stat_t_df"]] <- data.frame(subtype = levels(as.factor(tmp_df$subtype)))#one sample t-test
    for(i in 1:5){
      tmp_list[["perc_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["perc_df"]],1,
                                                         FUN = function(x){
                                                           sum(cell_meta$cell_type == x[["subtype"]] &
                                                                 cell_meta$interval == paste0(i,"0um"))/
                                                             sum(cell_meta$interval == paste0(i,"0um") &
                                                                   cell_meta$cell_type %in% cell_type_list[[j]])})
      tmp_list[["dens_df"]][[paste0("R",i*10)]]  <- apply(tmp_list[["dens_df"]],1,
                                                          FUN = function(x){
                                                            1000000 * sum(cell_meta$cell_type == x[["subtype"]] &
                                                                            cell_meta$interval == paste0(i,"0um"))/
                                                              (dilate_vec[[paste0(i,"0um")]]/res_factor^2) })
      
      tmp_list[["stat_t_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["stat_t_df"]],1,
                                                           FUN = function(x){
                                                             tmp_vec2 <- as.numeric(tmp_df[[paste0("R",i*10)]]) %>% 
                                                               .[tmp_df$subtype == x[["subtype"]]] %>% na.omit()
                                                             
                                                             tmp_obj <- t.test(x = tmp_vec2,
                                                                               mu = tmp_vec[[x[["subtype"]]]], 
                                                                               alternative = ifelse(mean(tmp_vec2) > tmp_vec[[x[["subtype"]]]],
                                                                                                    "greater", "less"))
                                                             tmp_obj$p.value
                                                           })
      
      tmp_list[["stat_c_df"]][[paste0("R",i*10)]] <- apply(tmp_list[["stat_c_df"]],1,
                                                           FUN = function(x){
                                                             a = sum(cell_meta$cell_type == x[["subtype"]] & cell_meta$interval == paste0(i,"0um"))
                                                             b = sum(cell_meta$cell_type != x[["subtype"]] & cell_meta$interval == paste0(i,"0um"))
                                                             c = sum(cell_meta$cell_type == x[["subtype"]] & cell_meta$interval != paste0(i,"0um"))
                                                             d = sum(cell_meta$cell_type != x[["subtype"]] & cell_meta$interval != paste0(i,"0um"))
                                                             
                                                             tmp_obj <- suppressWarnings(chisq.test(x = matrix(c(a,c,b,d), nrow = 2)))
                                                             tmp_obj$p.value
                                                           })
    }
    #Add overall 
    tmp_list[["dens_df"]] <- left_join(tmp_list[["dens_df"]],
                                       data.frame(subtype = names(tmp_vec),overall = tmp_vec), by = "subtype")
    tmp_list[["dens_df"]][["overall"]] <- tmp_vec
    tmp_list[["perc_df"]][["overall"]] <- apply(tmp_list[["perc_df"]],1,
                                                FUN = function(x){
                                                  sum(cell_meta$cell_type == x[["subtype"]])/
                                                    sum(cell_meta$cell_type %in% cell_type_list[[j]])})
    return_ls[["stat_ls"]][[j]] <- tmp_list
    #Stacked Barplot
    tmp_df %<>% melt(id.vars = c("plaque_id","subtype")) %>% subset(!is.na(value))
    #Skip for CADG
    for(i in levels(factor(tmp_df$subtype))){
      if(sum(tmp_df$value == 0 & tmp_df$subtype == i)/sum(tmp_df$subtype == i) < 0.99){
        tmp_df <- tmp_df[tmp_df$subtype != i|tmp_df$value < quantile(tmp_df$value[tmp_df$subtype == i], probs = 0.99),]
      }
    }
    
    return_ls[["plot_ls"]][[j]] <- stacked_vln(tmp_df, cell_type_palette_list[[j]], j, "subtype", tmp_vec)
  }
  return_ls
}
