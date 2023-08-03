#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Subsitute missing by -1
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
h5n <- args[1]
prefix <- args[2]
setPrefix <- gsub("assays.h5", "", h5n)

## Load libraries
library(HDF5Array)
library(SummarizedExperiment)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

impute_matrix <- function(mat){
  
  medians <- matrixStats::rowMedians(mat, na.rm = T)
  na.mat <- which(is.na(mat), arr.ind = T)
  mat[na.mat] <- medians[na.mat[, 1]]
  mat
}

mat <- data.matrix(assay(SE))

message("Imputing matrix")
imp_mat <- impute_matrix(mat)

assay(SE) <- DelayedArray(imp_mat)

saveHDF5SummarizedExperiment(SE, ".", prefix = paste0(prefix, "missingSub_"))
