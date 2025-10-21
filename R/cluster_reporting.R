library(party)      # for ctree
library(FactoMineR) # for catdes
library(factoextra) # for the fviz_* functions
library(tidyverse)
library(glue)
library(writexl)
library(ComplexHeatmap)

#' Utility Function:  Unlisting any nested structure (except data frames)
## ----------------------------------------------------------------------------------
soft_unlist <- function(x){
  if (!("list" %in% class(x))){
    list(x)
  } else {
    unlist(map(x, soft_unlist), recursive = F)
  }
}

#' Produce a table showing the size of every cluster, given a
#' tibble mapping subject IDs to clusters.
#'
#' @param clusters A tibble containing a column named "Cluster"
#'
#' @return A table giving the size of each cluster
#' @export
#'
#' @examples
cluster_table <- function(clusters){
  table(clusters$Cluster)
}

#' Create an MCA scatterplot, with dots colored by cluster membership.
#'
#' @param clusters A tibble containing a column named "Cluster"
#' @param mca The output of FactoMineR::MCA(), applied to the data
#' that have been clustered
#' @param axes The axes to be plotted
#'
#' @return A ggplot object showing the selected dimensions of the MCA,
#' with points colored according to cluster membership.
#' @export
#'
#' @examples
cluster_plot <- function(clusters, mca, axes = c(1, 2), title = ""){
  plot <- fviz_mca_ind(mca, axes = axes, geom = "point",
                       habillage = as.factor(clusters$Cluster),
                       addEllipses = T, title = title)
}

#' Do random forest analysis based on a clustering
#'
#' @param df The data that have been clustered
#' @param clusters A tibble containing a column named "Cluster"
#'
#' @return A classification tree output by party::ctree(), with
#' the fewest levels necessary to distinguish among all clusters.
#' @export
#'
#' @examples
rf_from_clustering <- function(df, clusters){
  maxdepth <- 1 + ceiling(log2(max(as.numeric(clusters$Cluster))))
  tcm_rf <- ctree(formula = Cluster ~ ., data = df %>%
                    left_join(clusters) %>%
                    mutate(Cluster = as.factor(Cluster)) %>%
                    select(-Subject)#,
                  # controls =
                  #   ctree_control(maxdepth = maxdepth)
                  )
  return(tcm_rf)
}

#' Write a classification tree to a PDF file.
#'
#' @param rf The output of rf_from_clustering() (or party::ctree())
#' @param outputFile The path to a PDF file to use as output
#' @param width Width of the output file, in inches
#' @param height Height of the output file, in inches
#'
#' @return No return value, but plots the classification tree in the
#' specified output file.
#' @export
#'
#' @examples
export_rf <- function(rf, outputFile, width=20, height=10){
  pdf(outputFile, width = width, height = height)
  plot(rf)
  dev.off()
}

#' Get category descriptions (using v tests) for a given clustering
#'
#' @param df The data that have been clustered
#' @param clusters A tibble containing a column named "Cluster"
#'
#' @return A list containing the results of v tests for all variables
#' that have a significant association (positive or negative)
#' with each cluster
#' @export
#'
#' @examples
get_catdes <- function(df, clusters, proba = 0.01){
  catdes <- df %>%
    select(-Subject) %>%
    mutate(Cluster = as.factor(clusters$Cluster)) %>%
    {catdes(., ncol(.), proba = proba)}
  return(catdes)
}
#' Write the output of catdes to an Excel file, and plot the heatmap
#' visualization in a PDF file
#'
#' @param catdes The output of get_catdes() (or FactoMineR::catdes())
#' @param clusters A tibble containing a column named "Cluster"
#' @param pdfFile The path to a PDF file for plotting the catdes
#' heatmap
#' @param excelFile The path to an Excel file for writing the catdes
#' tables
#'
#' @return No return value, but plots a heatmap visualization of the
#' catdes results to the specified PDF file, and writes the catdes
#' tables to the specified Excel file (one per sheet).
#' @export
#'
#' @examples
# export_catdes <- function(catdes, clusters, pdfFile, excelFile){
#   pdf(pdfFile, width = 10, height = 5)
#   plot(catdes, sort = as.character(which.max(cluster_table(clusters))))
#   dev.off()
# 
#   write_xlsx(catdes %>%
#                soft_unlist() %>%
#                map(~ as.data.frame(.) %>%
#                      rownames_to_column("Variable")),
#              excelFile)
# }
export_catdes <- function(catdes, clusters, pdfFile, excelFile){
  pdf(pdfFile, width = 10, height = 10)
  plot(catdes %>% .$category %>% 
    imap(~ as.data.frame(..1) %>% 
           rownames_to_column("variable") %>% 
           filter(!str_detect(variable, "FALSE")) %>% 
           select(variable, v.test) %>% 
           mutate(cluster = ..2)) %>% 
    reduce(rbind) %>% 
    pivot_wider(names_from = cluster, values_from = v.test, values_fill = 0) %>% 
    column_to_rownames("variable") %>% as.matrix() %>% t() %>% 
    Heatmap(show_row_dend = F, 
            show_column_dend = F, 
            row_names_side = "left",
            row_names_max_width = unit(Inf, "cm"),
            column_names_side = "top",
            column_names_max_height = unit(Inf, "cm"),
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20),
            color_space = "RGB",
            col = c("red", "white", "blue"),
            cluster_rows = F,
            name = "V-Test"))
  plot(catdes %>% .$category %>% 
         imap(~ as.data.frame(..1) %>% 
                rownames_to_column("variable") %>% 
                filter(!str_detect(variable, "FALSE|^Region|AgeOnset|Combined")) %>% 
                select(variable, v.test) %>% 
                mutate(cluster = ..2)) %>% 
         reduce(rbind) %>% 
         pivot_wider(names_from = cluster, values_from = v.test, values_fill = 0) %>% 
         column_to_rownames("variable") %>% as.matrix() %>% t() %>% 
         .[,apply(., 2, function(x){max(abs(x))}) > 10] %>% 
         Heatmap(show_row_dend = F, 
                 show_column_dend = F, 
                 row_names_side = "left",
                 row_names_max_width = unit(Inf, "cm"),
                 column_names_side = "top",
                 column_names_max_height = unit(Inf, "cm"),
                 row_names_gp = gpar(fontsize = 20),
                 column_names_gp = gpar(fontsize = 20),
                 color_space = "RGB",
                 col = c("red", "white", "blue"),
                 cluster_rows = F,
                 name = "V-Test"))
  dev.off()
  
  write_xlsx(catdes %>%
               soft_unlist() %>%
               map(~ as.data.frame(.) %>%
                     rownames_to_column("Variable")),
             excelFile)
}
