#' ---
#' title: "TICgen_BHC"
#' output: html_notebook
#' ---
#' 4/12/2023 - Analysis of TICgen data using Bayesian Hierarchical Clustering
#' put in a fix line# 43-45
#' 6/26/2023 - added tcnft
#' 7/11/2023 - Added tcnfctd - CTD + TD with no flags + Region
#' 7/24/2023 - Added ctdnf - CTD + TD with no flags
#' 8/5/2024  - Added tcnftu - TD, NoFlags, One/Family
#' 6/12/2025 - Added tcf - tcnftu w/o race or ethnicity
#' 8/26/2025 - Added TD - One TD Child / Family with no flags, no missing Dx and
#'              date of diagnosis as proxy for missing TD and OC Age of Onset.
## ----------------------------------------------------------------------------------
library(doParallel)
library(readxl)
library(writexl)
library(cluster)
library(party)
library(FactoMineR)
library(factoextra)
library(kmed)
library(tidyverse)
library(BHC)
library(glue)
library(ComplexHeatmap)

#' Read in preprocessed data
## ----------------------------------------------------------------------------------
load("./data/KMeans.rdata")


run_bhc <- function(df) {
  dfmat <- df %>%
    column_to_rownames("Subject") %>%
    modify(as.numeric) %>% as.matrix()
  if (is.null(rownames(dfmat))) {
    rownames(dfmat) <- pull(df,Subject)
  }
  res_BHC <- bhc(dfmat,
      rownames(dfmat),
      robust = 1, randomised = T, m = 200, verbose = T,
      numThreads = detectCores() - 1)
  return(res_BHC)
}

get_bhc_clusters <- function(bhc, outputFile, df) {
  WriteOutClusterLabels(bhc, outputFile, verbose = T)
  clusters <- read_tsv(outputFile, col_types = "c") %>%
    transmute(Subject = x,
              Cluster =
                map_chr(x, ~ str_match(., "CLUSTER ([0-9]+)")[2]) %>% 
                as.numeric()) %>%
    fill(Cluster) %>%
    filter(!str_detect(Subject, "^---CLUSTER")) %>% 
    left_join(df %>% select(Subject), .)
  return(clusters)
}

source("./R/cluster_reporting.R")

write_bhc_reports <- function(df, clusters, mca, bhc_res){
  name <- deparse(substitute(df))
  
  message(glue("Writing results for optimal BHC of {name}."))
  
  if (!dir.exists(name)) {dir.create(name)}
  
  write_tsv(clusters, glue("{name}/{name}_BHC_clusters.tsv"))
  
  sink(glue("./output/{name}_BHC/{name}_BHC_cluster_sizes.txt"))
  print(cluster_table(clusters))
  sink()
  
  pdf(glue("./output/{name}_BHC/{name}_BHC_cluster_plots.pdf"), 5, 5)
  cluster_plot(clusters, mca, 
               title = glue("BHC clusters for {name}")) %>% plot()
  cluster_plot(clusters, mca, 
               title = glue("BHC clusters for {name}"), 
               axes = c(1,3)) %>% plot()
  cluster_plot(clusters, mca, 
               title = glue("BHC clusters for {name}"), 
               axes = c(2,3)) %>% plot()
  dev.off()
  
  rf <- rf_from_clustering(df, clusters)
  export_rf(rf, glue("./output/{name}_BHC/{name}_BHC_random_forest.pdf"))
  
  catdes <- get_catdes(df, clusters)
  export_catdes(catdes, clusters,
                glue("./output/{name}_BHC/{name}_BHC_catdes.pdf"),
                glue("./output/{name}_BHC/{name}_BHC_catdes.xlsx"))
  
  pdf(glue("./output/{name}_BHC/{name}_BHC_res.pdf"), 200, 20)
  plot(bhc_res)
  dev.off()
}

### TD - TD, no flags, One/Family
TD_res_BHC <- run_bhc(TD)
TD_clusters <- get_bhc_clusters(TD_res_BHC, 
                                    "TICgen_TD_BHC_labels.txt",
                                    TD)
write_bhc_reports(TD, TD_clusters, TD_mca, TD_res_BHC)
#

save(list = ls(),
     file = "./data/bhc.RData")
