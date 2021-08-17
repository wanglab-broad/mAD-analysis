#Get Gene Expr in 0-40um filtered by pct
min_pct <- 0.05
starmap_obj <- starmap
starmap_obj@misc <- list()

#Combined Clustering
i = "AD_mouse9494"
j = "AD_mouse9723"

gene_list <- union(row.names(starmap_obj@assays[["RNA"]]@counts)[
  rowSums(starmap_obj@assays[["RNA"]]@counts[,starmap_obj@meta.data[["sample"]] == i] > 0)
  > sum(starmap_obj@meta.data[["sample"]] == i) * min_pct],
  row.names(starmap_obj@assays[["RNA"]]@counts)[
    rowSums(starmap_obj@assays[["RNA"]]@counts[,starmap_obj@meta.data[["sample"]] == j] > 0)
    > sum(starmap_obj@meta.data[["sample"]] == j) * min_pct]
  )

tmp_df <- matrix(0, nrow = 2*length(gene_list), ncol = 5) 
tmp_df <- data.frame(gene = gene_list, sample = c(rep(i,length(gene_list)),
                                rep(j,length(gene_list))),tmp_df)
colnames(tmp_df)[3:6] <- paste("R",1:4,"0",sep = "")
colnames(tmp_df)[7] <- paste("R","40+",sep = "")
#Raw Counts
#tmp_mat <- starmap_obj@assays[["RNA"]]@counts[gene_list,starmap_obj@meta.data[["sample"]] == sample_name]
#Scaled
tmp_mat <- cbind(starmap_obj@assays[["scaled"]]@data[gene_list,starmap_obj@meta.data[["sample"]] == i],
                 starmap_obj@assays[["scaled"]]@data[gene_list,starmap_obj@meta.data[["sample"]] == j])
cell_meta <- rbind(data.frame(sample = i,starmap@misc$AD_mouse9494_plaque_all$AD_mouse9494_cell_meta),
                   data.frame(sample = j,starmap@misc$AD_mouse9723_plaque_all$AD_mouse9723_cell_meta))
for (k in 1:4) {
  tmp_df[1:length(gene_list),k+2] <- rowSums(tmp_mat[,cell_meta$interval == paste(k,"0um",sep = "") & cell_meta$sample == i]) /
    sum(cell_meta$interval == paste(k,"0um",sep = "") & cell_meta$sample == i)
  tmp_df[(length(gene_list)+1):(2*length(gene_list)),k+2] <- rowSums(tmp_mat[,cell_meta$interval == paste(k,"0um",sep = "") & 
                                                                         cell_meta$sample == j]) /
    sum(cell_meta$interval == paste(k,"0um",sep = "") & cell_meta$sample == j)
}
tmp_df[1:length(gene_list),7] <- rowSums(tmp_mat[,cell_meta$interval %in% c("50um","other") & cell_meta$sample == i]) /
  sum(cell_meta$interval %in% c("50um","other") & cell_meta$sample == i)
tmp_df[(length(gene_list)+1):(2*length(gene_list)),7] <- rowSums(tmp_mat[,cell_meta$interval %in% c("50um","other") & cell_meta$sample == j]) / 
  sum(cell_meta$interval %in% c("50um","other") & cell_meta$sample == j)

tmp_df[,3:7] <- t(scale(t(as.matrix(tmp_df[,3:7]))))
#Weighted averaging
#for(k in 1:4){
#  tmp_df[1:length(gene_list),k+2] <- tmp_df[1:length(gene_list),k+2] * 
#    sum(cell_meta$interval == paste(k,"0um",sep = "") & cell_meta$sample == i) + 
#    tmp_df[(length(gene_list)+1):(2*length(gene_list)),k+2] * 
#    sum(cell_meta$interval == paste(k,"0um",sep = "") & cell_meta$sample == j)
#}
#tmp_df[1:length(gene_list),7] <- tmp_df[1:length(gene_list),7] * 
#  sum(cell_meta$interval == "other" & cell_meta$sample == i) + 
#  tmp_df[(length(gene_list)+1):(2*length(gene_list)),7] * 
#  sum(cell_meta$interval == "other" & cell_meta$sample == j)


#kmeans try
library(factoextra)
fviz_nbclust(tmp_df[,3:7],kmeans,k.max = 15,method = "wss")
set.seed(2021)
tmp_df.km_res <- kmeans(tmp_df[,3:7],6,nstart = 25)
#PC
#fviz_cluster(tmp_df.km_res,t(scale(t(as.matrix(tmp_df)))),axes = c(1,2))

plot_list <- list()
for(k in 1:6){
  plot_list[[i]] <- Heatmap(tmp_df[tmp_df.km_res$cluster == k, 3:7],
                            cluster_columns = F,cluster_rows  = T, show_row_names = F,
                            column_title = paste(sample_name,min_pct,"Cluster",k)) 
  print(which(colSums(tmp_df[tmp_df.km_res$cluster == k, 3:7]) == max(colSums(tmp_df[tmp_df.km_res$cluster == k, 3:7]))))
}

#heatmap_seq <- c(6,3,1,5,4,2,7)
heatmap_seq <- c(6,4,5,3,2,1)
#row.names(tmp_df) <- str_to_title(row.names(tmp_df))
v_heatmap <- NULL
sample_name <- "AD_mouse9494"
for(k in heatmap_seq){
  v_heatmap <- v_heatmap %v% Heatmap(tmp_df %>% subset(subset = sample == sample_name) %>%
                                       .[tmp_df.km_res$cluster[tmp_df$sample == sample_name] == k, 3:7],
                                     cluster_columns = F,cluster_rows  = T, show_row_names = F, show_row_dend = F,
                                     row_title = paste("C",k),border_gp = gpar(col = "black", lty = 1, lwd = 2),
                                     column_title = paste(sample_name,"PCT:",min_pct,"P-val",p_val)) 
}
draw(v_heatmap, ht_gap = unit(0.1,"cm"))
#Stat test
p_val <- 0.01
gene_df <- data.frame(tmp_df[,1:2], cluster = tmp_df.km_res$cluster, p_val = 1, p_adj = 1)


#13mo: 1-4/5 2-2 3-1 4-3 5-!1 6-!1
#8mo: 1-3/4 2-2 3-3 4-4/5 5-1
gene_df$p_val <- apply(gene_df,1,
                       FUN = function(x){
                         i <- switch(x[["cluster"]],
                                     "6" = 1, "4" = 2, "5" = 3, "3" = 4, "2" = 5, "1" = 1)
                         tmp_vec1 <- tmp_mat[x[["gene"]], cell_meta$interval == paste(i,"0um",sep = "") & 
                                               cell_meta$sample == x[["sample"]]]
                         tmp_vec2 <- tmp_mat[x[["gene"]], cell_meta$interval != paste(i,"0um",sep = "") & 
                                               cell_meta$sample == x[["sample"]]]
                         a <- wilcox.test(tmp_vec1,tmp_vec2, alternative = "greater")$p.value
                         return(a)
                       })
gene_df$p_adj <- p.adjust(gene_df$p_val)
#gene_df$sig_interval <- apply(gene_df,1,
#                       FUN = function(x){
#                         tmp_vec <- rep(1,3)
#                         for(i in 1:3){
#                           tmp_vec1 <- tmp_mat[x[["gene"]],cell_meta$interval == paste(i,"0um",sep = "")]
#                           tmp_vec2 <- tmp_mat[x[["gene"]],cell_meta$interval != paste(i,"0um",sep = "")]
#                           tmp_vec[i] <- wilcox.test(tmp_vec1,tmp_vec2, alternative = "greater")$p.value
#                         }
#                         if(min(tmp_vec) < 0.01) return(which(tmp_vec == min(tmp_vec)))
#                         else return(0)
#                       })

v_heatmap <- NULL
sample_name <- "AD_mouse9494"
for(k in heatmap_seq[1:3]){
  v_heatmap <- v_heatmap %v% Heatmap(tmp_df %>% subset(subset = sample == sample_name) %>% `rownames<-`(str_to_title(.[,1])) %>%
                                       .[tmp_df.km_res$cluster[tmp_df$sample == sample_name] == k & 
                                           gene_df$p_val[tmp_df$sample == sample_name] < p_val, 3:7],
                                     cluster_columns = F,cluster_rows  = T, show_row_names = T, show_row_dend = T,
                                     row_title = paste("C",k),rect_gp = gpar(col = "black", lwd = 1),
                                     column_title = paste(sample_name,"PCT:",min_pct,"P-val",p_val)) 
}
draw(v_heatmap, ht_gap = unit(0.1,"cm"))

#Save 13mo
gene_cl_res_13 <- list()
gene_cl_res_13[["p_val"]] <- p_val
gene_cl_res_13[["expr_df"]] <- tmp_df
gene_cl_res_13[["km_res"]] <- tmp_df.km_res
gene_cl_res_13[["stat_test"]] <- gene_df[1:(1*length(gene_list)),]
gene_cl_res_13[["combined_heatmap"]] <- v_heatmap
gene_cl_res_13[["filtered_heatmap"]] <- v_heatmap
#reload
tmp_df <- gene_cl_res_13[["expr_df"]]   
tmp_df.km_res <- gene_cl_res_13[["km_res"]] 
gene_df <- gene_cl_res_13[["stat_test"]] 
#save 8mo
gene_cl_res_8 <- list()
gene_cl_res_8[["p_val"]] <- p_val
gene_cl_res_8[["expr_df"]] <- tmp_df
gene_cl_res_8[["km_res"]] <- tmp_df.km_res
gene_cl_res_8[["stat_test"]] <- gene_df[(length(gene_list)+1):(2*length(gene_list)),]
gene_cl_res_8[["combined_heatmap"]] <- v_heatmap
gene_cl_res_8[["filtered_heatmap"]] <- v_heatmap

#Export Heatmaps
pdf(file = "Aug10_combined_gene_cluster.pdf",width = 4,height = 10,useDingbats = FALSE)
gene_cl_res_8[["combined_heatmap"]] %>% print()
gene_cl_res_8[["filtered_heatmap"]] %>% print()
gene_cl_res_13[["combined_heatmap"]] %>% print()
gene_cl_res_13[["filtered_heatmap"]] %>% print()
dev.off()
#Draw Venn Plot

require(ggVennDiagram)

pdf(file = "Aug10_combined_Venn_diagram.pdf",width = 4,height = 4,useDingbats = FALSE)
tmp_vec1 <- as.character(gene_cl_res_8[["stat_test"]]$gene[gene_cl_res_8[["stat_test"]]$p_val < p_val &
                                                             gene_cl_res_8[["stat_test"]]$cluster %in% heatmap_seq[1:3]])
tmp_vec2 <- as.character(gene_cl_res_13[["stat_test"]]$gene[gene_cl_res_13[["stat_test"]]$p_val < p_val &
                                                              gene_cl_res_13[["stat_test"]]$cluster %in% heatmap_seq[1:3]])
ggVennDiagram(list("8mo" = tmp_vec1,
                   "13mo" = tmp_vec2,
                   "PIG" = PIG))+labs(title = paste("Excluded 10µm Down","PCT:",min_pct,"P-val:",p_val)) %>% print()
dev.off()
#Export csv
gene_df$cluster <- as.factor(gene_df$cluster)
levels(gene_df$cluster) <- c("Cluster 6", "Cluster 5", "Cluster 4", "Cluster 2", "Cluster 3", "Cluster 1")
write_csv(data.frame(gene_df[,1:4],tmp_df[,3:7]) %>% 
            subset(subset = sample == "AD_mouse9494") %>% 
            mutate(gene = str_to_title(gene), cluster = s.character(cluster)) %>% 
            arrange(cluster) %>% .[,-2],
          path = "Aug11_gene_cl_13mo.csv")

#Export Heatmpa by clust
for(i in c(1:6)){
  Heatmap(t(scale(t(tmp_df[tmp_df.km_res$cluster == i & gene_df$p_val < 0.05, 1:5]))),
          cluster_columns = F,cluster_rows  = T, show_row_names = T,
          row_title = paste("C",i), row_names_gp = gpar(fontsize = 5),
          column_title = paste(sample_name,"PCT:",min_pct,"P-val",p_val)) %>% draw()
}
dev.off()

#PCA Visualization
tmp_df.pca <- prcomp(t(scale(t(as.matrix(gene_cl_res_13[["expr_df"]])))))
tmp_df.pca_val <- get_pca_ind(tmp_df.pca)

plotly::plot_ly(x = tmp_df.pca_val$coord[,1],
                y = tmp_df.pca_val$coord[,2],
                z = tmp_df.pca_val$coord[,3],
                type = "scatter3d",mode="markers",
                color = as.factor(gene_cl_res_13$km_res$cluster),
                size = 5,
                alpha = 0.8)

cl_colors <- brewer.pal(6,"Set1")
cl_colors <- cl_colors[as.numeric(gene_cl_res_13$km_res$cluster)]
pdf("gene_cluster_pca.pdf", width = 8,height = 6,useDingbats = FALSE)
s3d <- scatterplot3d(tmp_df.pca_val$coord[, 1:3], pch = 16, grid=T, box=T,
                     color = cl_colors,main = "8months")
legend(s3d$xyz.convert(1, 4, 1), legend = levels(as.factor(gene_cl_res_13$km_res$cluster)),
       col =  brewer.pal(6,"Set2"), pch = 16)

dev.off()

#Validation Heatmap
gene_list <- as.character(gene_cl_res_8$stat_test$gene[gene_cl_res_8$stat_test$cluster %in% c(2,3,5,6) & 
                                                         gene_cl_res_8$stat_test$p_val < 0.01])
gene_list <- as.character(gene_cl_res_13$stat_test$gene[gene_cl_res_13$stat_test$cluster %in% c(1,2,4,5) & 
                                                          gene_cl_res_13$stat_test$p_val < 0.01])
cell_meta <- starmap_64@misc$AD_mouse9721_plaque_all$AD_mouse9721_cell_meta
sample_name <- "AD_mouse9721"
cell_meta <- starmap_64@misc$AD_mouse9919_plaque_all$AD_mouse9919_cell_meta
sample_name <- "AD_mouse9919"

gene_list <- intersect(gene_list,row.names(starmap_64@assays[["RNA"]]@counts))

tmp_df <- matrix(0, nrow = length(gene_list), ncol = 5) %>% as.data.frame(row.names = gene_list)
colnames(tmp_df)[1:4] <- paste("R",1:4,"0",sep = "")
colnames(tmp_df)[5] <- paste("R","40+",sep = "")

tmp_mat <- starmap_64@assays[["RNA"]]@counts[gene_list,starmap_64@meta.data$sample == sample_name]
for (i in 1:4) {
  tmp_df[,i] <- rowSums(tmp_mat[,cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33]) /
    sum(cell_meta$min_border_dist < i * 33 & cell_meta$min_border_dist >= (i-1) * 33)
}
tmp_df[,5] <- rowSums(tmp_mat[,cell_meta$min_border_dist > 4 * 33]) / sum(cell_meta$min_border_dist > 4 * 33)

row.names(tmp_df) <- str_to_title(row.names(tmp_df))

Heatmap(t(scale(t(tmp_df))),
        cluster_columns = F,cluster_rows  = T, show_row_names = T,
        column_title = paste(sample_name,"Validation")) 

fviz_pca_ind(tmp_df.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
#Hierarchical clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

# Dissimilarity matrix
d <- dist(t(scale(t(as.matrix(tmp_df[,1:4])))), method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

sub_grp <- cutree(hc1, k = 4)

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

plot_list <- list()



Heatmap(t(scale(t(tmp_df[,1:4]))),
        cluster_columns = F,cluster_rows  = T, show_row_names = F,column_title = paste(sample_name,min_pct)) 


#PIG list

PIG <- c('APOE','AXL','CD63','CD63-PS','CD9','CTSB','CTSD','CTSL','CTSZ','H2-K1','HEXA','LGALS3BP','LYZ2',
         'NPC2','TREM2','TYROBP','H2-D1','B2M','C4B','GFAP','SERPINA3N','ARPC1B','C1QA','C1QB','C1QC','C4A',
         'CLU','CSF1R','CST3','CTSA','CTSS','CTSH','CX3CR1','CYBA','FCER1G','FCGR3','FCRLS','GRN','GUSB','GNS',
         'GPX4','HEXB','IGFBP5','ITGB5','ITM2B','LAPTM5','LGMN','LY86','MAN2B1','MPEG1','OLFML3','PLEK','PRDX6',
         'GPX4-PS','S100A6','RPL18A','VSIR')


