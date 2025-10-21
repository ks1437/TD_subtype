###
# KMeans.r
###
library(tidyverse)
library(glue)
library(parallel)

# Sourcing functions for finding optimal K means
source("./R/optimalK.r")

# Reporting using the functions in optimalK.R and cluster_reporting.R
source("./R/cluster_reporting.R")

###
# Functions
###
binarise <- function(df){
  return(df %>%
           {cbind(select(., Subject),
                  select_if(., is.logical),
                  select_if(select(., -Subject), is.factor) %>%
                    tab.disjonctif())} %>%
           modify_if(is.numeric, as.logical))
}

write_mca_report <- function(mca) {
  name <- deparse(substitute(mca))
  
  write_csv(as_tibble(mca$eig), file = glue("./output/{name}_eigenvalues.csv"))
  write_csv(as_tibble(mca$var$contrib), file = glue("./output/{name}_contrib.csv"))
  write_csv(as_tibble(mca$var$cos2), file = glue("./output/{name}_cos2.csv"))
}

scale_data <- function(df){
  df_scaled <- df %>%
    modify_if(is.factor, as.numeric) %>%
    column_to_rownames("Subject") %>%
    select_if(~ length(table(.)) > 1) %>%
    modify(scale)
  return(df_scaled)
}

get_optimal_clusters <- function(df, df_scaled){
  optimal_k_res <- find_optimal_k(df_scaled)
  clusters <- optimal_k_res$optimal_kmeans %>%
    {cbind(df %>% select(Subject),
           tibble(Cluster = .$cluster))}
  return(clusters)
}

reduce_mca_dims <- function(mca){
  name <- deparse(substitute(mca))
  
  optimal_ks <- get_kmeans_stability(mca, nrep = 100)
  
  pdf(glue("./output/{name}_stability_plot.pdf"), width = 10, height = 5)
  plot_kmeans_stability(optimal_ks) %>% plot()
  dev.off()
  
  ndim <- find_optimal_ndim(optimal_ks)
  
  mca_coords <- mca$ind$coord[, 1:ndim] %>%
    scale()
  
  return(mca_coords)
}

write_kmeans_reports <- function(df, coords, clusters, mca, name) {
  message(glue("Writing results for optimal K-means clustering of {name}."))
  
  if (!dir.exists(name)) {
    dir.create(paste0("./output/",name))
  }
  
  message("Writing elbow plot...")
  pdf(glue("./output/{name}/{name}_KMeans_elbow_plot.pdf"), 5, 5)
  elbow_plot(coords) %>% plot()
  dev.off()
  
  message("Writing silhouette plot...")
  pdf(glue("./output/{name}/{name}_KMeans_silhouette_plot.pdf"), 5, 5)
  silhouette_plot(coords) %>% plot()
  dev.off()
  
  message("Writing table of cluster sizes...")
  sink(glue("./output/{name}/{name}_KMeans_cluster_sizes.txt"))
  print(cluster_table(clusters))
  sink()
  
  message("Plotting MCA with points colored by cluster membership...")
  pdf(glue("./output/{name}/{name}_KMeans_cluster_plots.pdf"), 5, 5)
  cluster_plot(clusters, mca) %>% plot()
  cluster_plot(clusters, mca, axes = c(1,3)) %>% plot()
  cluster_plot(clusters, mca, axes = c(2,3)) %>% plot()
  dev.off()
  
  message("Writing output of random forest...")
  rf <- rf_from_clustering(df, clusters)
  export_rf(rf, glue("./output/{name}/{name}_KMeans_random_forest.pdf"))
  
  message("Writing output of catdes...")
  catdes <- get_catdes(df, clusters)
  export_catdes(catdes, clusters,
                glue("./output/{name}/{name}_KMeans_catdes.pdf"),
                glue("./output/{name}/{name}_KMeans_catdes.xlsx"))
}
###################################################

# load("./data/preproc.RData")
load ("./data/kmeans.RData")

###
# Process TD
###
### TD - Child with TD, one per family, no flags or unrateable Dx

TD_binarised <- binarise(TD)

TD_mca <- MCA(TD_binarised %>%
                 select(-Subject), ncp = Inf)

write_mca_report(TD_mca)

# %variance explained by n dimensions
cumsum(TD_mca$eig[,2])
ndim_kaiser <- max(which(TD_mca$eig[,1] > 0.1))
TD_mca_loadings <- TD_mca$svd$U
ndim_nontrivial <- colSums(abs(TD_mca$var$coord) > 0.3) %>%
  cummin() %>% `>`(3) %>% which() %>% max()
# Noise filter <= 97% variance
VAR_CUTOFF <- 70.0
ndim <- last(which(cumsum(TD_mca$eig[,2]) <= VAR_CUTOFF))

### Raw Data (Categorical) *** Not useful ***
TD_scaled <- TD %>%
  column_to_rownames("Subject") %>%
  modify_if(is.factor, as.numeric) %>%
  select_if(~ length(table(.)) > 1) %>%
  modify(scale)

TD_optimal_k_res <- find_optimal_k(TD_scaled)
TD_clusters <- TD_optimal_k_res$optimal_kmeans %>%
  {cbind(TD %>% select(Subject),
         tibble(Cluster = .$cluster))}
#############

TD_mca_optimal_ks <- get_kmeans_stability(TD_mca, nrep = 100)

TD_mca_ndim <- find_optimal_ndim(TD_mca_optimal_ks)

TD_mca_coords <- TD_mca$ind$coord[, 1:TD_mca_ndim] %>%
  scale()

###
# Applying K-means to (scaled) TD
write_kmeans_reports(TD, TD_scaled, TD_clusters, TD_mca, "TD")
pdf("./output/TD_mca/TD_mca_stability_plot.pdf", width = 10, height = 5)
plot_kmeans_stability(TD_mca_optimal_ks) %>% plot()
dev.off()
# Applying K-means to the MCA coordinates for TD, 5 clusters (optimal)
TD_mca5_clusters <-
  cbind(TD %>% select(Subject),
        tibble(Cluster = kmeans(TD_mca_coords,
                                centers = 5,
                                iter.max = ITER_MAX,
                                nstart = NSTART)$cluster))
write_kmeans_reports(TD, TD_mca_coords, TD_mca5_clusters, TD_mca, "TD_mca5")
###
# Applying K-means to (scaled) TD k=6
TD_mca6_clusters <-
  cbind(TD %>% select(Subject),
        tibble(Cluster = kmeans(TD_mca_coords,
                                centers = 6,
                                iter.max = ITER_MAX,
                                nstart = NSTART)$cluster))
# Applying K-means to the MCA coordinates for TD, 6 clusters
write_kmeans_reports(TD, TD_mca_coords, TD_mca6_clusters, TD_mca, "TD_mca6")
###
# Applying K-means to (scaled) TD k=7
TD_mca7_clusters <-
  cbind(TD %>% select(Subject),
        tibble(Cluster = kmeans(TD_mca_coords,
                                centers = 7,
                                iter.max = ITER_MAX,
                                nstart = NSTART)$cluster))
# Applying K-means to the MCA coordinates for TD, 7 clusters
write_kmeans_reports(TD, TD_mca_coords, TD_mca7_clusters, TD_mca, "TD_mca7")
###
# Adding vertical bar at cluster 7
###
x <- elbow_plot(TD_mca_coords)
x + geom_vline(xintercept=7,linetype='dashed', color='blue',size=0.5)
###
# end TD
###


save(list = ls(),
     file = "./data/kmeans.RData", version = 3)
#-------------------------------------------------------------#


###
# Try with 8 dimensions instead of 3 identified by cluster stability
###
# TD_mca_coords <- TD_mca$ind$coord[, 1:8] %>%
#   scale()
# ###
# TD_mca_optimal_k_res <- find_optimal_k(TD_mca_coords)
# TD_mca_clusters <- TD_mca_optimal_k_res$optimal_kmeans %>%
#   {cbind(TD %>% select(Subject),
#          tibble(Cluster = .$cluster))}
###
