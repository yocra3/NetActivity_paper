#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Output projects labels
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

project <- SE$project
if ("sample_type" %in% colnames(colData(SE))){
  project[SE$sample_type == "Solid Tissue Normal"] <- "Normal"
}

write.table(project, file = paste0(setPrefix, "_individuals_labels.txt"), quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

