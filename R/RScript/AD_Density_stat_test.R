#
library(imager)
sample_name <- "AD_mouse9494"
cell_type <- "Top_level"
region <- "Cortex"

top_level_palette_list <- c("#5fd9c3", "#ebeb00", "#29e322", "#b8adff", "#f2402c", "#5bb0eb", 
                            "#698270", "#f9cfff", "#fcb900", "#8600d4", "#c800cc", "#ab4b66","#b3b125","#8c8c8c")
names(top_level_palette_list) <- c("Ex", "Inhi", "Astro", "CA1", "CA2", "CA3", "DG", "Endo", "Micro", "OPC", "Oligo",
                                   "SMC","LHb","Mix")

i = 5
#cell_type_list and top_level_palette_list pre-defined
ad_density_stat_test <- function(starmap, sample_name, region){
  region_img <- starmap@misc[[paste(sample_name, "morph", sep = "_")]]$region
  cell_meta <- starmap@misc[[paste(sample_name,"plaque",region, sep = "_")]][[paste(sample_name,"cell_meta", sep = "_")]]
  cell_meta$cell_type[cell_meta$cell_type == "OPC_Gpr17"] <- "OPC"

  message("Generate Seed Image")
  plaque_meta <- starmap@misc[[paste(sample_name,"plaque",region, sep = "_")]][[paste(sample_name,"plaque_meta", sep = "_")]]
  plaque_img <- starmap@misc[[paste(sample_name,"plaque",region, sep = "_")]][[paste(sample_name,"plaque_img", sep = "_")]]
  plaque_seed_img <- matrix(0, nrow = dim(plaque_img)[1], ncol = dim(plaque_img)[2])
  for(j in 1:dim(plaque_meta)[1]){
    plaque_seed_img[round(plaque_meta$m.cx[j]),round(plaque_meta$m.cy[j])] <- plaque_meta$plaque_id[j]
  }
  plaque_seed_img <- EBImage::Image(plaque_seed_img)
  res_list <- list()
  for(j in c(names(cell_type_list),"Top_level")){
    if(j == "Top_level"){
      res_list[[j]] <- list(res_pval_df = data.frame(cell_type = names(top_level_palette_list),
                                                     row.names = names(top_level_palette_list)),
                            res_diff_df = data.frame(cell_type = names(top_level_palette_list),
                                                     row.names = names(top_level_palette_list)),
                            res_dens_df = data.frame(cell_type = names(top_level_palette_list),
                                                     row.names = names(top_level_palette_list)),
                            res_sdev_df = data.frame(cell_type = names(top_level_palette_list),
                                                     row.names = names(top_level_palette_list)),
                            res_cont_df = data.frame(cell_type = names(top_level_palette_list),
                                                     row.names = names(top_level_palette_list)))
    }else{
      res_list[[j]] <- list(res_pval_df = data.frame(cell_type = cell_type_list[[j]],row.names = cell_type_list[[j]]),
                            res_diff_df = data.frame(cell_type = cell_type_list[[j]],row.names = cell_type_list[[j]]),
                            res_dens_df = data.frame(cell_type = cell_type_list[[j]],row.names = cell_type_list[[j]]),
                            res_sdev_df = data.frame(cell_type = cell_type_list[[j]],row.names = cell_type_list[[j]]),
                            res_cont_df = data.frame(cell_type = cell_type_list[[j]],row.names = cell_type_list[[j]]))
    }
  }
  
  require(EBImage)
  plaque_dilate_img <- plaque_img
  plaque_dilate_img_bk <- plaque_img
  
  for(i in 1:5){
    message(paste("Dilate",i,"of",5,sep = " "))
    plaque_dilate_img_bk <- plaque_dilate_img
    plaque_dilate_img <- EBImage::dilate(plaque_img, kern = makeBrush(round(3.175*20*i), shape = "disc")) %>% as.matrix()
    
    if(region != "all"){
      x <- switch(region,
                  "Cortex" = "1",
                  "White Matter" = "2",
                  "Hippocampus" = "3")
      plaque_dilate_img <- plaque_dilate_img * (region_img == x)
    }
    
    
    plaque_ws_img <- propagate(plaque_dilate_img, plaque_seed_img,mask = plaque_dilate_img)
    
    if(i != 1) plaque_ws_img <- plaque_ws_img * (1-plaque_dilate_img_bk)
    
    
    
    plaque_dilate_meta <- computeFeatures.shape(plaque_ws_img) %>% 
      as.data.frame() %>% data.frame(., plaque_id = as.numeric(row.names(.)))
    plaque_dilate_meta_bk <- plaque_dilate_meta
    
    for(cell_type in c(names(cell_type_list),"Top_level")){
      plaque_dilate_meta <- plaque_dilate_meta_bk
      if(cell_type == "Top_level"){
        for(j in names(top_level_palette_list)){
          plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]] <- 0
        }
        level = "top_level"
      }else{
        for(j in cell_type_list[[cell_type]]){
          plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]] <- 0
        }
        level = "cell_type"
      }
      
      
      for(j in 1:dim(cell_meta)[1]){
        if(plaque_ws_img[cell_meta$x[j],cell_meta$y[j]] %in% plaque_dilate_meta$plaque_id &
           (cell_meta[[level]][j] %in% names(top_level_palette_list) | cell_meta[[level]][j] %in% cell_type_list[[cell_type]])){
          #print(paste(as.character(cell_meta$cell_type[j]),j))
          tmp <- as.numeric(plaque_ws_img[cell_meta$x[j],cell_meta$y[j]])
          plaque_dilate_meta[[paste(cell_meta[[level]][j],"_",i,"0um_count",sep = "")]][row.names(plaque_dilate_meta) == tmp] <- 1 +
            plaque_dilate_meta[[paste(cell_meta[[level]][j],"_",i,"0um_count",sep = "")]][row.names(plaque_dilate_meta) == tmp] 
        }
      }
      if(cell_type == "Top_level"){
        for(j in names(top_level_palette_list)){
          plaque_dilate_meta[[paste(j,"_",i,"0um_density",sep = "")]] <- 
            1000000 * plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]]/plaque_dilate_meta$s.area*3.175^2
        }
        density_vec <- rep(0,length(names(top_level_palette_list))); names(density_vec) <- names(top_level_palette_list)
      }else{
        for(j in cell_type_list[[cell_type]]){
          plaque_dilate_meta[[paste(j,"_",i,"0um_density",sep = "")]] <- 
            1000000 * plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]]/plaque_dilate_meta$s.area*3.175^2
        }
        density_vec <- rep(0,length(cell_type_list[[cell_type]])); names(density_vec) <- cell_type_list[[cell_type]]
      }
      
      
      
      res_pval_df <- res_list[[cell_type]][["res_pval_df"]]
      res_diff_df <- res_list[[cell_type]][["res_diff_df"]]
      res_dens_df <- res_list[[cell_type]][["res_dens_df"]]
      res_sdev_df <- res_list[[cell_type]][["res_sdev_df"]]
      res_cont_df <- res_list[[cell_type]][["res_cont_df"]]
      res_pval_df[[paste(i,"0um",sep = "")]] <- 0
      res_diff_df[[paste(i,"0um",sep = "")]] <- 0
      res_dens_df[[paste(i,"0um",sep = "")]] <- 0
      res_sdev_df[[paste(i,"0um",sep = "")]] <- 0
      res_cont_df[[paste(i,"0um",sep = "")]] <- 0
      
      for(j in names(density_vec)){
        if(region != "all"){
          x <- switch(region,
                      "Cortex" = sum(region_img == "1"),
                      "White Matter" = sum(region_img == "2"),
                      "Hippocampus" = sum(region_img == "3"))
        }else x <- dim(plaque_img)[1] * dim(plaque_img)[2]

        density_vec[[j]] <- 1000000 * sum(cell_meta[[level]] == j)/x*3.175^2
        res_sdev_df[j,paste(i,"0um",sep = "")] <- sd(plaque_dilate_meta[[paste(j,"_",i,"0um_density",sep = "")]])
        res_dens_df[j,paste(i,"0um",sep = "")] <- 1000000 * (sum(plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]])/
                                                               sum(plaque_dilate_meta[["s.area"]])*3.175^2)
        res_cont_df[j,paste(i,"0um",sep = "")] <- sum(plaque_dilate_meta[[paste(j,"_",i,"0um_count",sep = "")]])
        res_diff_df[j,paste(i,"0um",sep = "")] <- res_dens_df[j,paste(i,"0um",sep = "")] - density_vec[[j]]
        #res_pval_df[j,paste(i,"0um",sep = "")] <- wilcox.test(x = plaque_dilate_meta[[paste(j,"_",i,"0um_density",sep = "")]],
        #                                                      mu = density_vec[[j]],
        #                                                      alternative = ifelse(res_diff_df[j,paste(i,"0um",sep = "")] > 0,
        #                                                                       "greater","less"))$p.value
        res_pval_df[j,paste(i,"0um",sep = "")] <- t.test(x = plaque_dilate_meta[[paste(j,"_",i,"0um_density",sep = "")]],
                                                              mu = density_vec[[j]], 
                                                              alternative = ifelse(res_diff_df[j,paste(i,"0um",sep = "")] > 0,
                                                                               "greater","less"))$p.value
        
      }
      if(i == 5) res_dens_df[["overall"]] <- density_vec
      res_list[[cell_type]][[paste(i,"0um_dilate_df",sep = "")]] <- plaque_dilate_meta
      res_list[[cell_type]][["res_pval_df"]] <- res_pval_df
      res_list[[cell_type]][["res_diff_df"]] <- res_diff_df
      res_list[[cell_type]][["res_dens_df"]] <- res_dens_df
      res_list[[cell_type]][["res_sdev_df"]] <- res_sdev_df
      res_list[[cell_type]][["res_cont_df"]] <- res_cont_df
    }
 
  }
  return(res_list)
}
tmp <- list()
tmp[["AD9494_C"]] <- ad_density_stat_test(starmap, "AD_mouse9494", "Cortex")
tmp[["AD9494_H"]] <- ad_density_stat_test(starmap, "AD_mouse9494", "Hippocampus")
tmp[["AD9723_C"]] <- ad_density_stat_test(starmap, "AD_mouse9723", "Cortex")
tmp[["AD9723_H"]] <- ad_density_stat_test(starmap, "AD_mouse9723", "Hippocampus")
tmp[["AD9494_A"]] <- ad_density_stat_test(starmap, "AD_mouse9494", "all")
tmp[["AD9723_A"]] <- ad_density_stat_test(starmap, "AD_mouse9723", "all")
tmp[["AD9919_A"]] <- ad_density_stat_test(starmap_64, "AD_mouse9919", "all")
tmp[["AD9721_A"]] <- ad_density_stat_test(starmap_64, "AD_mouse9721", "all")
res_list <- tmp[["AD9723_C"]]

res_pval_df <- res_list[["Ex"]][["res_pval_df"]]
res_sdev_df <- res_list[["Ex"]][["res_sdev_df"]]

for(i in res_list){
  res_pval_df <- rbind(res_pval_df,i[["res_pval_df"]])
  res_sdev_df <- rbind(res_sdev_df,i[["res_sdev_df"]])
}
res_pval_df <- res_pval_df[-c(5:8),]
res_sdev_df <- res_sdev_df[-c(5:8),]

write_csv(res_sdev_df,"Jul31_AD9723_Cortex_density_sdev.csv")
write_csv(res_pval_df,"Jul31_AD9723_Cortex_density_pval_wilcox.csv")

ad_plot_cell_dens <- function(density_df,
                               palette_list, #named list contains color and name of each cluster
                               title = ""){
  tmp_df <- density_df %>% melt(id.vars = "cell_type")
  tmp_df$cell_type <- factor(tmp_df$cell_type,levels = names(palette_list))
  g <- ggplot(data = tmp_df,aes(x = variable, y = value, fill = .data[["cell_type"]]))+
    geom_bar(stat="identity") + 
    labs(title = title) +
    scale_fill_manual(values = palette_list,
                      breaks = names(palette_list),
                      labels = names(palette_list)) #+ 
    #theme_classic()
  
  
  return(g)
}

#Re-plot sub-type density
plot_list <- list()
#
for(i in c("Micro","Astro","Oligo","Inhi","Top_level")){
  for(j in c(9494,9723)){
    for(k in c("C","H")){
      plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
        ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]],
                          palette_list = cell_type_palette_list[[i]],
                          title = paste("AD",j,k,i,sep = "_"))
    }
  }
}
i = "Ex"
k = "C"
for(j in c(9494,9723)){
  plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
    ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]],
                      palette_list = cell_type_palette_list[[i]],
                      title = paste("AD",j,k,i,sep = "_"))
}
i = "CADG"
k = "H"
for(j in c(9494,9723)){
  plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
    ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]],
                      palette_list = cell_type_palette_list[[i]],
                      title = paste("AD",j,k,i,sep = "_"))
}
i = "Top_level"
for(j in c(9494,9723)){
  for(k in c("C","H")){
    plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
      ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]] %>% subset(!cell_type %in% c("DG","Mix","LHb")),
                        palette_list = top_level_palette_list[-c(7,13,14)],#2766
                        title = paste("AD",j,k,i,sep = "_"))
  }
}
for(j in c(9494,9723)){
  k = "A"
  plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
    ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]] %>% subset(!cell_type %in% c("DG","Mix")),
                      palette_list = top_level_palette_list[-c(7,14)],#2766
                      title = paste("AD",j,k,i,sep = "_"))
}
for(j in c(9919,9721)){
  for(k in c("A")){
    plot_list[[paste("AD",j,k,i,sep = "_")]] <- 
      ad_plot_cell_dens(tmp[[paste("AD",j,"_",k,sep = "")]][[i]][["res_dens_df"]] %>% subset(!cell_type %in% c("DG","LHb")),
                        palette_list = top_level_palette_list[-c(7,13)],#64
                        title = paste("AD",j,k,i,sep = "_"))
  }
}
pdf(file = "Aug04_subtype_density_theme.pdf",width = 12,height = 9,useDingbats = FALSE)
plot_list %>% print()
dev.off()





pdf(file = "Aug02_64gene_toplevel_density.pdf",width = 12,height = 9,useDingbats = FALSE)
plot_list %>% print()
dev.off()
#Export CSV
cell_stat_res <- tmp
#Plaque Stat
for(i in c("AD9494_C","AD9494_H","AD9723_C","AD9723_H")){
  for(j in names(cell_type_list)){
    for(k in cell_type_list[[j]]){
      tmp_df <- cell_stat_res[[i]][[j]][["10um_dilate_df"]][,1:7]
      for(a in 1:5){
        require(dplyr)
        tmp_df <- left_join(tmp_df,
                            cell_stat_res[[i]][[j]][[paste(a,"0um_dilate_df",sep = "")]][,c("plaque_id",
                                                                                            paste(k,"_",a,"0um_count",sep = ""))],
                            by = "plaque_id")
          
      }
      for(a in 1:5){
        require(dplyr)
        tmp_df <- left_join(tmp_df,
                            cell_stat_res[[i]][[j]][[paste(a,"0um_dilate_df",sep = "")]][,c("plaque_id",
                                                                                            paste(k,"_",a,"0um_density",sep = ""))],
                            by = "plaque_id")
        
      }
      write_csv(tmp_df,path = paste("stat_res/",i,"/",i,"_",j,"_",str_replace(k,"/","-"),".csv",sep = ""))
    }
  }
}


#Dens+wilcox
for(i in c("Ex","Inhi","Micro","Astro","Oligo","CADG")){
  tmp_df <- data.frame(region = rep(c(rep("Cortex",length(cell_type_list[[i]])),
                                        rep("Hippocampus",length(cell_type_list[[i]]))),4),
                       rbind(cell_stat_res[["AD9723_C"]][[i]][["res_dens_df"]],cell_stat_res[["AD9723_H"]][[i]][["res_dens_df"]],
                             cell_stat_res[["AD9494_C"]][[i]][["res_dens_df"]],cell_stat_res[["AD9494_H"]][[i]][["res_dens_df"]],
                             cbind(rbind(cell_stat_res[["AD9723_C"]][[i]][["res_pval_df"]],
                                              cell_stat_res[["AD9723_H"]][[i]][["res_pval_df"]],
                                              cell_stat_res[["AD9494_C"]][[i]][["res_pval_df"]],
                                              cell_stat_res[["AD9494_H"]][[i]][["res_pval_df"]]),
                                   data.frame(overall = rep("-",length(cell_type_list[[i]])*4)))))
  write_csv(tmp_df,path = paste("stat_res/",i,"_dens_n_ttest.csv",sep = ""))
}
#Top Level
tmp_df <- data.frame(region = rep(c(rep("All",11),rep("Cortex",11),
                                    rep("Hippocampus",11)),4),
                     rbind(cell_stat_res[["AD9723_A"]][[i]][["res_dens_df"]][tmp_vec,],
                           cell_stat_res[["AD9723_C"]][[i]][["res_dens_df"]][tmp_vec,],
                           cell_stat_res[["AD9723_H"]][[i]][["res_dens_df"]][tmp_vec,],
                           cell_stat_res[["AD9494_A"]][[i]][["res_dens_df"]][tmp_vec,],
                           cell_stat_res[["AD9494_C"]][[i]][["res_dens_df"]][tmp_vec,],
                           cell_stat_res[["AD9494_H"]][[i]][["res_dens_df"]][tmp_vec,],
                           cbind(rbind(cell_stat_res[["AD9723_A"]][[i]][["res_pval_df"]][tmp_vec,],
                                       cell_stat_res[["AD9723_C"]][[i]][["res_pval_df"]][tmp_vec,],
                                       cell_stat_res[["AD9723_H"]][[i]][["res_pval_df"]][tmp_vec,],
                                       cell_stat_res[["AD9494_A"]][[i]][["res_pval_df"]][tmp_vec,],
                                       cell_stat_res[["AD9494_C"]][[i]][["res_pval_df"]][tmp_vec,],
                                       cell_stat_res[["AD9494_H"]][[i]][["res_pval_df"]][tmp_vec,]),
                                 data.frame(overall = rep("-",11*6)))))
write_csv(tmp_df,path = "stat_res/Top_level_dens_n_ttest.csv")


#Percentage+Chi-square Test
for(i in c("Ex","Inhi","Micro","Astro","Oligo","CADG")){
  tmp_df <- data.frame(
    sample = c(rep("AD_mouse9723",length(cell_type_list[[i]])*2),rep("AD_mouse9494",length(cell_type_list[[i]])*2)),
    region = rep(c(rep("Cortex",length(cell_type_list[[i]])),rep("Hippocampus",length(cell_type_list[[i]]))),2),
    rbind(cell_stat_res[["AD9723_C"]][[i]][["res_cont_df"]],
          cell_stat_res[["AD9723_H"]][[i]][["res_cont_df"]],
          cell_stat_res[["AD9494_C"]][[i]][["res_cont_df"]],
          cell_stat_res[["AD9494_H"]][[i]][["res_cont_df"]]),
    overall = 0)
                       
  tmp_df[["overall"]] <- apply(tmp_df,1,
                               FUN = function(x){
                                 return(sum(starmap@meta.data$cell_type == x[["cell_type"]] & 
                                            starmap@meta.data$sample == x[["sample"]] &
                                            starmap@meta.data$region == x[["region"]]))
                               })
  tmp_df2 <- tmp_df
  for(j in 1:5){
    for(k in 1:dim(tmp_df)[1]){
      a <- switch(paste(tmp_df2$sample[k],tmp_df2$region[k],sep = "_"),
                  "AD_mouse9723_Cortex" = 
                    sum(cell_stat_res[["AD9723_C"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                  "AD_mouse9723_Hippocampus" = 
                    sum(cell_stat_res[["AD9723_H"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                  "AD_mouse9494_Cortex" = 
                    sum(cell_stat_res[["AD9494_C"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                  "AD_mouse9494_Hippocampus" = 
                    sum(cell_stat_res[["AD9494_H"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]))
      a <- a - tmp_df[k,j+3]
      tmp_df2[k,j+3] <- chisq.test(matrix(c(tmp_df[k,j+3],a,tmp_df[k,9] - tmp_df[k,j+3],
                                          sum(starmap@meta.data$sample == as.character(tmp_df2$sample[k]) &
                                                starmap@meta.data$region == as.character(tmp_df2$region[k])) -
                                            a - tmp_df[k,9]),ncol = 2))$p.value
    }
  }
  for(j in c("AD_mouse9494","AD_mouse9723")){
    for(k in c("Cortex","Hippocampus")){
      tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9] <- tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9]/
        matrix(rep(colSums(tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9]),length(cell_type_list[[i]])), 
               nrow = length(cell_type_list[[i]]), byrow = T) * 100
    }
  }
  write_csv(rbind(tmp_df,tmp_df2),path = paste("stat_res/",i,"_perc_n_chisq.csv",sep = ""))
}

#Top level
tmp_df <- data.frame(
  sample = c(rep("AD_mouse9723",length(tmp_vec)*3),rep("AD_mouse9494",length(tmp_vec)*3)),
  region = rep(c(rep("All",length(tmp_vec)),rep("Cortex",length(tmp_vec)),rep("Hippocampus",length(tmp_vec))),2),
  rbind(cell_stat_res[["AD9723_A"]][[i]][["res_cont_df"]][tmp_vec,],
        cell_stat_res[["AD9723_C"]][[i]][["res_cont_df"]][tmp_vec,],
        cell_stat_res[["AD9723_H"]][[i]][["res_cont_df"]][tmp_vec,],
        cell_stat_res[["AD9494_A"]][[i]][["res_cont_df"]][tmp_vec,],
        cell_stat_res[["AD9494_C"]][[i]][["res_cont_df"]][tmp_vec,],
        cell_stat_res[["AD9494_H"]][[i]][["res_cont_df"]][tmp_vec,]),
  overall = 0)

tmp_df[["overall"]] <- apply(tmp_df,1,
                             FUN = function(x){
                               if(x[["region"]] == "All"){
                                 return(sum(starmap@meta.data$top_level == x[["cell_type"]] & 
                                              starmap@meta.data$sample == x[["sample"]]))
                               }else{
                                 return(sum(starmap@meta.data$top_level == x[["cell_type"]] & 
                                              starmap@meta.data$sample == x[["sample"]] &
                                              starmap@meta.data$region == x[["region"]]))
                               }
                             })
tmp_df2 <- tmp_df
for(j in 1:5){
  for(k in 1:dim(tmp_df)[1]){
    a <- switch(paste(tmp_df2$sample[k],tmp_df2$region[k],sep = "_"),
                "AD_mouse9723_Cortex" = 
                  sum(cell_stat_res[["AD9723_C"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                "AD_mouse9723_Hippocampus" = 
                  sum(cell_stat_res[["AD9723_H"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                "AD_mouse9494_Cortex" = 
                  sum(cell_stat_res[["AD9494_C"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                "AD_mouse9494_Hippocampus" = 
                  sum(cell_stat_res[["AD9494_H"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                "AD_mouse9723_All" = 
                  sum(cell_stat_res[["AD9723_A"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]),
                "AD_mouse9494_All" = 
                  sum(cell_stat_res[["AD9494_A"]][["Top_level"]][["res_cont_df"]][[paste(j,"0um",sep = "")]]))
    a <- a - tmp_df[k,j+3]
    if(tmp_df2$region[k] == "All"){
      tmp_df2[k,j+3] <- chisq.test(matrix(c(tmp_df[k,j+3],a,tmp_df[k,9] - tmp_df[k,j+3],
                                            sum(starmap@meta.data$sample == as.character(tmp_df2$sample[k])) -
                                              a - tmp_df[k,9]),ncol = 2))$p.value
    }else{
      tmp_df2[k,j+3] <- chisq.test(matrix(c(tmp_df[k,j+3],a,tmp_df[k,9] - tmp_df[k,j+3],
                                            sum(starmap@meta.data$sample == as.character(tmp_df2$sample[k]) &
                                                  starmap@meta.data$region == as.character(tmp_df2$region[k])) -
                                              a - tmp_df[k,9]),ncol = 2))$p.value
    }
    
  }
}
for(j in c("AD_mouse9494","AD_mouse9723")){
  for(k in c("Cortex","Hippocampus","All")){
    tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9] <- tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9]/
      matrix(rep(colSums(tmp_df[tmp_df$sample == j & tmp_df$region == k,4:9]),length(tmp_vec)), 
             nrow = length(tmp_vec), byrow = T) * 100
  }
}
write_csv(rbind(tmp_df,tmp_df2),path = paste("stat_res/",i,"_perc_n_chisq.csv",sep = ""))



