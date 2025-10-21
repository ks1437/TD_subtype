###
# Functions for determining optimal number of clusters
# optimalK.r - New methods from STDHA books
#' 8/15/2024 - updated cluster stability, handles nulls
#' 1/2/2025 - Fixed problem of doParallel in nested loops 
###
library(cluster)
library(NbClust)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(foreach)
library(furrr)
library(party)
library(glue)

options(mc.cores = parallel::detectCores() - 1)

# load("./data/preproc.RData")

set.seed(15423)

# nboot = 50 to keep the function speedy.
# recommended value: nboot= 500 for your analysis.
KMAX <- 20
NBOOT <- 500
NSTART <- 1000 # Make this much larger
ITER_MAX <- 5000
# Noise filter <= 97% variance
VAR_CUTOFF <- 97.0
NREP <- 100

# tcm_mca_coords <- tcm_mca$ind$coord[, 1:last(which(cumsum(tcm_mca$eig[,2]) < VAR_CUTOFF))]
# #  tcm_mca$ind$coord[, 1:5]
# df <- tcm_mca_coords %>%
#   scale()

# Elbow method - K-Means - Optimal K
# Doesn't automatically tell you where the elbow is.
elbow_plot <- function(df){
  fviz_nbclust(df, kmeans, method = "wss", k.max = KMAX) +
    labs(subtitle = "Elbow method")
}

# Silhouette method
silhouette_plot <- function(df){
  fviz_nbclust(df, kmeans, method = "silhouette", k.max = KMAX,
               iter.max = ITER_MAX, nstart = NSTART)+
    labs(subtitle = "Silhouette method")
}

#' Finding optimal k for k-means.
## ----------------------------------------------------------------------------------------
find_optimal_k <- function(df, max_k = KMAX,
                           iter.max = ITER_MAX, nstart = NSTART){
  df_scaled <- df %>%
    # column_to_rownames("Subject") %>%
    # select_if(~ length(unique(.)) > 1) %>%
    # tab.disjonctif.prop() %>%
    scale()

  # plan(multisession, workers = detectCores() - 1)
  all_kmeans <- future_map(2:max_k, function(k){
    df_scaled %>%
      kmeans(., centers = k, iter.max = iter.max, nstart = nstart)
  }, .options = furrr_options(seed = TRUE))
  silhouettes <- future_map_dbl(all_kmeans, function(x){df_scaled %>%
      {silhouette(x$cluster, dist(df_scaled))} %>% summary() %>% .$avg.width})
  # plan(sequential)

  optimal_k <- which.max(silhouettes) + 1
  return(list(optimal_k = optimal_k,
              optimal_kmeans = all_kmeans[[optimal_k - 1]],
              silhouettes = silhouettes))
}

# Finding optimal number of MCA dimensions for stable K-means
# get_kmeans_stability <- function(mca, nrep = 10, max_k = KMAX,
get_kmeans_stability <- function(mca, nrep = 100, max_k = KMAX,
                                 iter.max = ITER_MAX, nstart = NSTART){
  optimal_ks <- tibble(Dimensions = numeric(0),
                       `Optimal K` = numeric(0))
  furrr_options(seed = T)
  # plan(multisession, workers = availableCores(omit = 1))
  plan(multisession, workers = 22)
  for (i in 1:sum(mca$eig[,3] < 95)){
  # for (i in 1:sum(mca$eig[,3] < 70)){
    message(paste("Trying", i, "dimensions", collapse = " "))
    mca_subset <- mca$ind$coord[, 1:i] %>%
      scale()
    subset_optimal_ks <-
      future_map_dbl(1:nrep,
                     possibly(~ find_optimal_k(mca_subset,
                                               max_k = KMAX,
                                               iter.max = ITER_MAX,
                                               nstart = NSTART)$optimal_k,
                              otherwise = NA_integer_))
    optimal_ks <- optimal_ks %>%
      rbind(tibble(Dimensions = i, `Optimal K` = subset_optimal_ks))
  }
  plan(sequential)

  return(optimal_ks)
}

plot_kmeans_stability <- function(optimal_ks) {
  name <- deparse(substitute(optimal_ks)) %>% 
    str_remove_all("_optimal_ks$")
  ggplot(optimal_ks, aes(Dimensions, `Optimal K`)) +
    geom_count() +
    ggtitle(glue("Stability plot for {name}"))
}

find_optimal_ndim <- function(optimal_ks) {
  ndim <- optimal_ks %>% 
    filter(!is.na(`Optimal K`)) %>%
    distinct() %>%
    group_by(Dimensions) %>%
    count() %>%
    ungroup() %>%
    mutate(n = if_else(cumall(n > 1), 1L, n),
           stable = cumall(n == 1)) %>%
    filter(stable) %>%
    .$Dimensions %>% last()
  return(ndim)
}

# Gap statistic
# Use verbose = FALSE to hide computing progression.

# fviz_nbclust(df, kmeans, nstart = NSTART,
#              method = "gap_stat", nboot = NBOOT,
#              k.max = KMAX)+
#   labs(subtitle = "Gap statistic method")
# tcm_clusgap <- clusGap(df, kmeans, KMAX, NBOOT, nstart = NSTART)
# fviz_gap_stat(tcm_clusgap)

# K-Means using optimal_k
# tcm_clusters <- kmeans(df, optimal_k, iter.max = ITER_MAX, nstart = 1000) %>%
#   {tibble(Cluster = .$cluster) %>%
#       cbind(tcm %>% select(Subject))} %>%
#   select(Subject, Cluster)
# table(tcm_clusters$Cluster)

# ggplot(as.data.frame(df) %>%
#          mutate(Cluster = as.factor(tcm_clusters$cluster)),
#        aes(`Dim 1`, `Dim 2`, col = Cluster)) +
#   geom_point() +
#   coord_fixed()

# fviz_mca_ind(tcm_mca, axes = c(1, 2), geom = "point",
#              habillage = as.factor(tcm_clusters$cluster),
#              addEllipses = T)
#
# fviz_mca_ind(tcm_mca, axes = c(1, 3), geom = "point",
#              habillage = as.factor(tcm_clusters$cluster),
#              addEllipses = T)
#
# fviz_mca_ind(tcm_mca, axes = c(2, 3), geom = "point",
#              habillage = as.factor(tcm_clusters$cluster),
#              addEllipses = T)
#
# ggsave("tcm_mca_27dims_6clusters.pdf", width = 10, height = 10)

# tcm_catdes <- catdes(tcm %>%
# #                       .[complete.cases(.),] %>%
#                        select(-Subject) %>%
#                        mutate(Cluster = as.factor(tcm_clusters$cluster)),
#                      ncol(tcm))
# pdf("tcm_catdes.pdf", width = 10, height = 5)
# plot(tcm_catdes, level = 0.05,
#      sort = 1, cex.names = 0.5)
# dev.off()
#
# tcm_rf <- ctree(formula = Cluster ~ .,
#                 data = tcm %>%
#                   .[complete.cases(.),] %>%
#                   select(-Subject) %>%
#                   mutate(Cluster = as.factor(tcm_clusters$cluster)),
#                 controls =
#                   ctree_control(maxdepth = 1 + ceiling(log2(max(tcm_clusters$cluster)))))

###
# Nbclust - 30 methods to find optimal k
###
# MINK <- 2
# MAXK <- 20
# DIST <- c("euclidean","manhattan","maximum","canberra", "binary", "minkowski")
# CMETHOD <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
#              "median", "centroid", "kmeans")
#
# plan(multisession, workers = parallel::detectCores() - 1)
#
# all_fviz <- expand.grid(DIST, CMETHOD) %>%
#   rename(distance = Var1, cmethod = Var2) %>%
#   {future_map2(.$distance, .$cmethod,
#                possibly(~ NbClust(df, distance = .x, min.nc = MINK,
#                                   max.nc = MAXK, method = .y) %>%
#                           fviz_nbclust(), otherwise = NA))}
###