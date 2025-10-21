###
#' cluster_comparison.R - Compare KMeans and BHC clusters
#' Heatmap showing % overlap of membership between KMeans and BHC clusters
#' works best.  Rest are experimental.
#' 6/30/2024 - V1.0
#' 8/15/2024 - V1.1 - changed tcnft -> TD & tcnft_mca -> TD_mca6
#' 2/04/2025 - V1.2 - Get info from Rdata files and use rownorm
#' 10/7/2025 - V1.3 - Reformat row and col labels and change name.
###

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(factoextra)

load("./data/KMeans.RData")
load("./data/BHC.RData")

# TD_bhc_clusters <- read_xlsx("BHC/TD/TD_BHC_catdes.xlsx", sheet = "call.X")
# TD_kmeans_clusters <- read_xlsx("Kmeans/TD_mca6/TD_mca6_KMeans_catdes.xlsx", sheet = "call.X")

TD_clusters_combined <- 
  inner_join(TD_clusters %>% rename(BHC = Cluster),
             TD_mca5_clusters %>% rename(KMeans = Cluster)) %>% 
  select(-Subject)

# TD_clusters_combined <- tibble(BHC = TD_bhc_clusters$Cluster,
#                                   KMeans = TD_kmeans_clusters$Cluster)

# Labeled cross-tab heatmap  *** This works best ***
pdf(file = "./output/KM5_BHC_Cluster_Overlap_Heatmap_colnorm.pdf")
plot(bhc_kmeans_crosstab %>% 
       `rownames<-`(paste0(rownames(.), " (", rowSums(.), ")")) %>% 
       `colnames<-`(paste0(colnames(.), " (", colSums(.), ")")) %>% 
       sweep(2, colSums(.), "/") %>% 
       round(2) %>% sweep(2, colSums(.), "/") %>% 
       Heatmap(col = c("white", "red"), cell_fun = function(j, i, x, y, width, height, fill) {
         grid.rect(x = x, y = y, width = width, height = height, 
                   gp = gpar(col = "black", fill = NA))
         if (.[i, j] > 0){
           grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 20, fontface = "bold"))
         }
       }, show_row_dend = F, show_column_dend = F,
       row_title = "BHC clusters", column_title = "K-means clusters",
       row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
       row_order = order(bhc_kmeans_crosstab %>% sweep(2, colSums(.), "/") %>% apply(1, max), decreasing = T),
       column_order = order(bhc_kmeans_crosstab %>% sweep(2, colSums(.), "/") %>% apply(2, max), decreasing = T),
       column_names_rot = 0, column_names_centered = T,
       name = "overlap"))
dev.off()

pdf(file = "./output/KM5_BHC_Cluster_Overlap_Heatmap_rownorm.pdf")
plot(bhc_kmeans_crosstab %>% 
       `rownames<-`(paste0(rownames(.), " (", rowSums(.), ")")) %>% 
       `colnames<-`(paste0(colnames(.), " (", colSums(.), ")")) %>% 
       sweep(1, rowSums(.), "/") %>% 
       round(2) %>% sweep(1, rowSums(.), "/") %>% 
       Heatmap(col = c("white", "red"), cell_fun = function(j, i, x, y, width, height, fill) {
         grid.rect(x = x, y = y, width = width, height = height, 
                   gp = gpar(col = "black", fill = NA))
         if (.[i, j] > 0){
           grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 20, fontface = "bold"))
         }
       }, show_row_dend = F, show_column_dend = F,
       row_title = "BHC clusters", column_title = "K-means clusters",
       row_title_gp = gpar(fontsize = 20), column_title_gp = gpar(fontsize = 20),
       row_order = order(bhc_kmeans_crosstab %>% sweep(1, rowSums(.), "/") %>% apply(1, max), decreasing = T),
       column_order = order(bhc_kmeans_crosstab %>% sweep(1, rowSums(.), "/") %>% apply(2, max), decreasing = T),
       column_names_rot = 0, column_names_centered = T,
       name = "overlap"))
dev.off()
############
# combn(1:3, 2) %>% t() %>% as.data.frame() %>% 
#   {map2(.$V1, .$V2, function(x, y){fviz_mca_ind(TD_mca6, axes = c(x, y), geom = "point",
#                                  habillage = as.factor(TD_bhc_clusters$Cluster),
#                                  addEllipses = T, ellipse.type = "convex") + 
#       ggpubr::stat_chull(data = TD_mca6$ind$coord[,c(x, y)] %>% 
#                            `colnames<-`(c("x", "y")) %>% 
#                            as.data.frame() %>% 
#                            mutate(cluster = TD_kmeans_clusters$Cluster), 
#                          mapping = aes(group = cluster), linewidth = 1,
#                          alpha = 0.1, geom = "polygon")})}
# 
# combn(1:3, 2) %>% t() %>% as.data.frame() %>% 
#   {map2(.$V1, .$V2, ~ fviz_mca_ind(TD_mca6, axes = c(..1, ..2), geom = "point",
#                                    habillage = as.factor(TD_kmeans_clusters$Cluster),
#                                    addEllipses = T, ellipse.type = "convex"))}
# 
# fviz_mca_ind(TD_mca6, axes = c(1, 2), geom = "point",
#              habillage = as.factor(TD_bhc_clusters$Cluster),
#              addEllipses = T, ellipse.type = "convex") + 
#   stat_ellipse(data = TD_mca6$ind$coord[,1:2] %>% 
#                  `colnames<-`(c("x", "y")) %>% 
#                  as.data.frame() %>% 
#                  mutate(cluster = TD_kmeans_clusters$Cluster), 
#                mapping = aes(group = cluster))
# 
# ggplot(as.data.frame(TD_mca6$ind$coord) %>% 
#          mutate(cluster_bhc = as.factor(TD_bhc_clusters$Cluster)) %>% 
#          group_by(cluster_bhc) %>% 
#          mutate(x = mean(`Dim 1`), y = mean(`Dim 2`)),
#        mapping = aes(`Dim 1`, `Dim 2`, col = cluster_bhc, fill = cluster_bhc)) + 
#   geom_point() + 
#   ggpubr::stat_chull(geom = "polygon", alpha = 0.1) + 
#   geom_label(aes(x = x, y = y, label = cluster_bhc), fill = "white", size = 8) + 
#   coord_fixed() + 
#   theme_minimal()
# 
# ggplot(as.data.frame(TD_mca6$ind$coord) %>% 
#          mutate(cluster_kmeans = TD_kmeans_clusters$Cluster,
#                 cluster_bhc = TD_bhc_clusters$Cluster) %>% 
#          group_by(cluster_kmeans, cluster_bhc) %>% 
#          mutate(n = n()) %>% 
#          group_by(cluster_kmeans) %>% 
#          mutate(cluster_bhc = if_else(n/n() > 0.1, cluster_bhc, NA_integer_)) %>% 
#          ungroup() %>% 
#          modify_if(is.integer, as.factor),
#        mapping = aes(`Dim 1`, `Dim 2`, col = cluster_kmeans, fill = cluster_kmeans)) + 
#   geom_point() + 
#   ggpubr::stat_chull(geom = "polygon", alpha = 0.1) + 
#   ggpubr::stat_chull(mapping = aes(group = cluster_bhc), geom = "polygon", col = "black", alpha = 0.1) + 
#   facet_wrap(vars(cluster_kmeans)) + 
#   coord_fixed() + 
#   theme_minimal()

### Doesn't look good.  Not used
# bhc_kmeans_crosstab <- table(TD_clusters_combined)
# Heatmap(bhc_kmeans_crosstab, col = c("white", "red"))
###