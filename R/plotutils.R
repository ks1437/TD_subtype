###
# plotutils.R - Create visualizations
# 08/21/2025 - V1.0 - Data Preprocessing and Exploration
###

library(rcompanion)
library(ComplexHeatmap)
library(flexmix)
library(glue)
library(tidyverse)
###
# Functions
###
pairwise_chisq <- function(df) {
  chisq <- outer(colnames(df), colnames(df), Vectorize(function(x, y) {
    chisq.test(table(df[[x]], df[[y]]))$statistic
  })) %>%
    `rownames<-`(colnames(df)) %>%
    `colnames<-`(colnames(df))
  cramerv <- outer(colnames(df), colnames(df), Vectorize(function(x, y) {
    cramerV(table(df[[x]], df[[y]]), bias.correct = T)
  })) %>%
    `rownames<-`(colnames(df)) %>%
    `colnames<-`(colnames(df))
  pvals <- outer(colnames(df), colnames(df), Vectorize(function(x, y) {
    chisq.test(table(df[[x]], df[[y]]))$p.value
  })) %>%
    `rownames<-`(colnames(df)) %>%
    `colnames<-`(colnames(df))
  pvals[lower.tri(pvals)] <- NA
  pvals <- p.adjust(pvals) %>%
    matrix(ncol(df), ncol(df), dimnames = list(colnames(df), colnames(df)))
  pvals[lower.tri(pvals)] <- t(pvals)[lower.tri(pvals)]
  return(list(
    statistic = chisq,
    cramerv = cramerv,
    p.value = pvals
  ))
}
###
triangularize <- function(x) {
  x[upper.tri(x, T)] <- NA
  return(x)
}
###
############################## End Functions #############################

plotParamCorr <- function(df) {
  dsname <- deparse(substitute(df))
  print("##########################################################")
  print(paste0("Cramer's V Param Correlation for: ", dsname))
  print("##########################################################")
  
  df_chisq <- pairwise_chisq(df %>% select(-Subject))$statistic
  df_chisq.pvals <- pairwise_chisq(df %>% select(-Subject))$p.value
  df_cramer <- pairwise_chisq(df %>% select(-Subject))$cramer

  pdf(file = paste0("./output/",dsname,"_ParamCorr_CramersV.pdf"))  
  Heatmap(
    df_cramer,
    show_row_dend = F,
    show_column_dend = F,
    col = c("white", "red"),
    column_title = paste0("Cramer's V for ",dsname, " columns"),
    name = "V-Test",
    # The following argument puts the Cramer's V values in the cells
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        gp = gpar(col = "black", fill = NA)
      )
      if (df_cramer[i, j] > 0) {
        grid.text(sprintf("%.1f", df_cramer[i, j]), x, y, gp = gpar(fontsize = 10))
      }
    }
  )
  dev.off()
  
  df_chisq_long <- inner_join(
    df_chisq %>%
      triangularize() %>%
      as.data.frame() %>%
      rownames_to_column("Var1") %>%
      pivot_longer(-Var1, names_to = "Var2", values_to = "X2") %>%
      arrange(desc(abs(X2))) %>%
      filter(!is.na(X2)),
    df_chisq.pvals %>%
      triangularize() %>%
      as.data.frame() %>%
      rownames_to_column("Var1") %>%
      pivot_longer(-Var1, names_to = "Var2", values_to = "p.value")
  ) %>%
    filter(p.value <= 0.05)
  
  write_csv(df_chisq_long, paste0("./output/",dsname,"_chisq.csv"))
}