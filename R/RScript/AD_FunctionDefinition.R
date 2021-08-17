#####Required Packages#####
library(Seurat)
library(SeuratDisk)
require(raster)
require(stringr)
require(dplyr)
require(parallel)
require(pbapply)
require(ggplot2)
library(readr)
require(EBImage)
#####Initialize#####
ad_ini <- function(){
  require(Seurat)
  require(SeuratDisk)
  require(raster)
  require(stringr)
  require(dplyr)
  require(parallel)
  require(pbapply)
  require(ggplot2)
  require(readr)
  require(EBImage)
  require(ggpubr)
  require(monocle3)
  print("Initialization Finished")
}


#####Cell Stat#####
#Test Parameter
#cell_meta_all <- starmap@misc[["AD_mouse9494_plaque_all"]][["AD_mouse9494_cell_meta"]]
#cell_meta_c   <- starmap@misc[["AD_mouse9494_plaque_Cortex"]][["AD_mouse9494_cell_meta"]]
#cell_meta_wm  <- starmap@misc[["AD_mouse9494_plaque_White Matter"]][["AD_mouse9494_cell_meta"]]
#cell_meta_h   <- starmap@misc[["AD_mouse9494_plaque_Hippocampus"]][["AD_mouse9494_cell_meta"]]
#
#img_size_all  <- sum(starmap@misc[["AD_mouse9494_morph"]][["region"]]  > 0)
#img_size_c    <- sum(starmap@misc[["AD_mouse9494_morph"]][["region"]] == 1)
#img_size_wm   <- sum(starmap@misc[["AD_mouse9494_morph"]][["region"]] == 2)
#img_size_h    <- sum(starmap@misc[["AD_mouse9494_morph"]][["region"]] == 3)
##
#res_df <- ad_cell_stat(starmap,"AD_mouse9494","top_level") #%>% write.csv(file = "AD9494_stat_toplevel.csv")
res_df <- ad_cell_stat(starmap,"AD_mouse9494","cell_type") #%>% write.csv(file = "AD9494_stat_sublevel.csv")
#res_df <- ad_cell_stat(starmap,"AD_mouse9723","top_level") #%>% write.csv(file = "AD9723_stat_toplevel.csv")
#res_df <- ad_cell_stat(starmap,"AD_mouse9723","cell_type") #%>% write.csv(file = "AD9723_stat_sublevel.csv")
#
#
##Export Latex Table
#res_df <- subset(res_df,top_level %in% c("Astro","Micro","Oligo","OPC") & region != "All")
#res_df$cell_type <- factor(as.character(res_df$cell_type),
#       levels = c('Astro', 'Astro_Cst3', 'Astro_Gfap/Vim', 'Micro', 'Micro_Gpr34', 'Micro_Cst7/Ctsb', 
#                  'Oligo', 'Oligo_Klk6','Oligo_Cldn11',  'OPC'),
#       labels = c('Astro1','Astro2','Astro3','Micro1','Micro2','Micro3','Oligo1','Oligo2','Oligo3','OPC'))
#res_df <- res_df %>% arrange(region)
#setwd("~/Desktop/At_Broad/AD/Export_csv/density_tables")
##Top 3 regions
#res_df <- subset(res_df, region == "All")
#res_df <- res_df %>% arrange(region)
#res_df[,c(1:2,3:8)] %>% write_csv(path = "AD9721_top_level_all_density.csv")
#res_df[,c(1:2,9:14)] %>% write_csv(path = "AD9723_top_level_all_percentage.csv")
##Cell type 
##Micro
#res_df2 <- subset(res_df, top_level == "Micro" & region != "All")
#res_df2[,c(1:3,4:9)] %>% write_csv(path = "AD9723_Micro_density.csv")
#res_df2 <- res_df2 %>% arrange(region)
#res_df2[,c(1:3,10:15)] %>% write_csv(path = "AD9723_Micro_percentage.csv")
##Astro
#res_df2 <- subset(res_df, top_level == "Astro" & region != "All")
#res_df2[,c(1:3,4:9)] %>% write_csv(path = "AD9723_Astro_density.csv")
#res_df2 <- res_df2 %>% arrange(region)
#res_df2[,c(1:3,10:15)] %>% write_csv(path = "AD9723_Astro_percentage.csv")
##Oligo
#res_df2 <- subset(res_df, top_level == "Oligo" & region != "All")
#res_df2[,c(1:3,4:9)] %>% write_csv(path = "AD9723_Oligo_density.csv")
#res_df2 <- res_df2 %>% arrange(region)
#res_df2[,c(1:3,10:15)] %>% write_csv(path = "AD9723_Oligo_percentage.csv")
##Ex
#res_df2 <- subset(res_df, top_level == "Ex" & region %in% c("Cortex","Hippocampus"))
#res_df2[,c(1:3,4:9)] %>% write_csv(path = "AD9723_Ex_density.csv")
#res_df2 <- res_df2 %>% arrange(region)
#res_df2[,c(1:3,10:15)] %>% write_csv(path = "AD9723_Ex_percentage.csv")
##Inhi
#res_df2 <- subset(res_df, top_level == "Inhi" & region %in% c("Cortex","Hippocampus"))
#res_df2[,c(1:3,4:9)] %>% write_csv(path = "AD9723_Inhi_density.csv")
#res_df2 <- res_df2 %>% arrange(region)
#res_df2[,c(1:3,10:15)] %>% write_csv(path = "AD9723_Inhi_percentage.csv")


ad_cell_stat <- function(starmap_obj,sample_name,level = "cell_type", only_all = F, raw_count = F){
  #rm Lhb SMC and combine OPC
  tmp_cell_meta_filter <- function(cell_meta){
    cell_meta <- subset(cell_meta, !top_level %in% c("SMC","LHb"))
    cell_meta[["cell_type"]][cell_meta[["cell_type"]] == "OPC_Gpr17"] <- "OPC"
    cell_meta
  }
  region_vec <- c("all","Cortex","White Matter","Hippocampus")
  region_alia_vec <- c("All","Cortex","Corpus callosum","Hippocampus")
  
  if(only_all){
    region_vec <- c("all")
    region_alia_vec <- c("All")
  }
  
  cell_meta_list <- list()
  for(i in region_vec){
    cell_meta_list[[i]] <- tmp_cell_meta_filter(starmap_obj@misc[[paste(sample_name,"plaque",i,sep = "_")]][[paste(sample_name,"cell_meta",sep = "_")]])
  }
  
  img_size_list <- c(0,0,0,0)
  names(img_size_list) <- region_vec
  for(i in region_vec){
    tmp <- switch(i,
                  "all" = 0,"Cortex"=1,"White Matter"=2,"Hippocampus"=3)
    if(i == "all") img_size_list[[i]] <- sum(starmap_obj@misc[[paste(sample_name,"morph",sep = "_")]][["plaque"]] >= 0)
    else img_size_list[[i]] <- sum(starmap_obj@misc[[paste(sample_name,"morph",sep = "_")]][["region"]] == tmp)
  }

  if(level == "cell_type"){
    tmp_df <- unique(cell_meta_list[["all"]][,c("top_level","cell_type")]) %>% arrange(top_level)
    res_df <- data.frame(top_level = rep(tmp_df[["top_level"]],each = length(region_vec)),
                         cell_type = rep(tmp_df[["cell_type"]],each = length(region_vec)),
                         #cell_type = rep(tmp_df[["cell_type_label"]],each = 4),
                         region = rep(region_alia_vec, dim(tmp_df)[1] ),
                         Density_10um = 0, Density_20um = 0, Density_30um = 0, Density_40um = 0, Density_50um = 0, Density_all = 0,
                         Percentage_10um = 0, Percentage_20um = 0, Percentage_30um = 0, Percentage_40um = 0, Percentage_50um = 0, Percentage_all = 0,
                         Pvalue_10um = 0,Pvalue_20um = 0,Pvalue_30um = 0,Pvalue_40um = 0,Pvalue_50um = 0)
  }else if(level == "top_level"){
    tmp_df <- unique(cell_meta_list[["all"]][,c("top_level")]) %>% sort()
    res_df <- data.frame(top_level = rep(tmp_df,each = length(region_vec)),
                         region = rep(region_alia_vec, length(tmp_df) ),
                         Density_10um = 0, Density_20um = 0, Density_30um = 0, Density_40um = 0, Density_50um = 0, Density_all = 0,
                         Percentage_10um = 0, Percentage_20um = 0, Percentage_30um = 0, Percentage_40um = 0, Percentage_50um = 0, Percentage_all = 0,
                         Pvalue_10um = 0,Pvalue_20um = 0,Pvalue_30um = 0,Pvalue_40um = 0,Pvalue_50um = 0)
  }
  
  tmp_cell_stat <- function(x,i,raw = F){
    region_name <- switch(x[["region"]],
                          "All" = "all","Cortex" = "Cortex","Corpus callosum" = "White Matter", "Hippocampus" = "Hippocampus")
    cell_meta <- cell_meta_list[[region_name]]
    if(i == "all"){
      region_size <- img_size_list[[region_name]]
    }else{
      region_size <- starmap_obj@misc[[paste(sample_name,"plaque",region_name,sep = "_")]][[paste(sample_name,"plaque_dilate_area",sep = "_")]][i] 
    }
     
    if(raw) region_size = 1
    else region_size = region_size/(3.3^2)
    if(i == "all"){
      return(sum(cell_meta[[level]] == x[[level]])/(region_size) )
    }else{
      return(sum(cell_meta[[level]] == x[[level]] & cell_meta[["min_border_dist"]] >= 33*(i-1) & cell_meta[["min_border_dist"]] < 33*i
                   )/(region_size))
    }
    
  }
  tmp_normalize <- function(vec){vec/sum(vec)}
  res_df$Density_all <- apply(res_df,1,
                              FUN = function(x){tmp_cell_stat(x,"all") * 1000000})
  res_df$Percentage_all <- apply(res_df,1,
                              FUN = function(x){tmp_cell_stat(x,"all",T)})
  

  if(level == "top_level" & !raw_count){
    for(i in region_alia_vec){
      res_df$Percentage_all[res_df$region == i] <- 100 * tmp_normalize(res_df$Percentage_all[res_df$region == i])
    }
  }else if(level == "cell_type" & !raw_count){
    cell_type_vec <- unique(cell_meta_list[["all"]][,c("top_level")]) %>% sort()
    for(k in cell_type_vec){
      for(i in region_alia_vec){
        res_df$Percentage_all[res_df$region == i & res_df$top_level == k] <- 
          100 * tmp_normalize(res_df$Percentage_all[res_df$region == i & res_df$top_level == k])
      }
    }
  }

  for(i in c(1:5)){
    res_df[[paste("Density_",i,"0um",sep = "")]] <- apply(res_df,1,
                                                          FUN = function(x){tmp_cell_stat(x,i) * 1000000})
    res_df[[paste("Percentage_",i,"0um",sep = "")]] <- apply(res_df,1,
                                                          FUN = function(x){tmp_cell_stat(x,i,T)})
    if(level == "top_level" & !raw_count){
      for(j in region_alia_vec){
        res_df[[paste("Percentage_",i,"0um",sep = "")]][res_df$region == j] <- 
          100 *tmp_normalize(res_df[[paste("Percentage_",i,"0um",sep = "")]][res_df$region == j])
      }
    }else if(level == "cell_type" & !raw_count){
      cell_type_vec <- unique(cell_meta_list[["all"]][,c("top_level")]) %>% sort()
      for(k in cell_type_vec){
        for(j in region_alia_vec){
          res_df[[paste("Percentage_",i,"0um",sep = "")]][res_df$region == j & res_df$top_level == k] <- 
            100 *tmp_normalize(res_df[[paste("Percentage_",i,"0um",sep = "")]][res_df$region == j & res_df$top_level == k])
        }
      }
    }
  }
  for(i in 1:5){
    res_df[[paste("Pvalue_",i,"0um",sep = "")]] <- apply(res_df,1,
                                                          FUN = function(x){
                                                            region_name <- switch(x[["region"]],
                                                                                  "All" = "all","Cortex" = "Cortex",
                                                                                  "Corpus callosum" = "White Matter", "Hippocampus" = "Hippocampus")
                                                            cell_meta <- cell_meta_list[[region_name]]
                                                            cell_partition_i    <- sum(cell_meta[[level]] == x[[level]] & 
                                                                                         cell_meta$min_border_dist >= (i-1)*33 &
                                                                                         cell_meta$min_border_dist <  i*33)
                                                            cell_partition_all  <- sum(cell_meta[[level]] == x[[level]]) - cell_partition_i
                                                            other_partition_i   <- sum(cell_meta[[level]] != x[[level]] & 
                                                                                         cell_meta$min_border_dist >= (i-1)*33 &
                                                                                         cell_meta$min_border_dist <  i*33)
                                                            other_partition_all <- sum(cell_meta[[level]] != x[[level]]) - other_partition_i
                                                            return(chisq.test(matrix(c(cell_partition_i, other_partition_i,
                                                                                       cell_partition_all, other_partition_all),
                                                                                     ncol = 2))$p.value)
                                                          })
  }
  return(res_df)
}
  
  
  
  
  
  
  
  
  
  
  






