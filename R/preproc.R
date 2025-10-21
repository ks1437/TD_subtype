###
#' preproc.R - Explore data and get break downs / splits
#' Date: 08/21/2025 - Initial version
###

library(pander)
library(readxl)
library(tidyverse)

SITE_REGION_XREF = "./input/SiteRegionList.xlsx"
TICGEN_CLIN_DATA = "./input/TICGen_Data_2023-07-15.xlsx"

# Get all utility functions
source("./R/stutils.R")
source("./R/plotutils.R")

tcg <- getInputData(SITE_REGION_XREF, TICGEN_CLIN_DATA)
tcx <- combineDiseaseDiag(tcg)
sink("./output/tcx_counts.txt")
reportSplits(tcx)
sink()
tct <- getOneChildPerFamWithTD(tcx)
sink("./output/tct_counts.txt")
reportDiagOnsetMultiSplits(tct)
sink()
tcno <- getNoOverrides(tcx)
sink("./output/tcno_counts.txt")
reportDiagOnsetMultiSplits(tcno)
sink()
tcnof <- getOneChildPerFamWithTDnoFlags(tcno)
sink("./output/tcnof_counts.txt")
reportDiagOnsetMultiSplits(tcnof)
sink()
# Remove columns that can't be used for clustering and rows (unratable)
TD <- genDataForClustering(tcnof)

###
# To Do: Sanity checks like # of diag = # of Age of Onset
###

write_tsv(TD, "./data/TD_Data.tsv")

save(list = ls(),
     file = "./data/preproc.RData", version = 3)
#-------------------------------------------------------------#

###############################################################
# TD Data Characteristics / Stats
###############################################################
###
# Write out TD dataset
###
write_tsv(TD, "./output/TD_Data.tsv")
###
# TD Stats
sink("./output/TD_counts.txt")
cat("Number of unique values in each column\n")
x <- TD %>% map_dbl( ~ length(table(.))) %>%
  as.data.frame()
x <- rownames_to_column(x)
colnames(x) <- c("Parameter", "Levels")
x %>% pander()
# Breakdown by category for each parameter
for (i in 2:ncol(TD)) {
  print (table(TD[i]) %>%
           as.data.frame() %>% pander())
}
sink()
###
# Plot Cramer's V Param Correlation
plotParamCorr(TD)
###
#####

