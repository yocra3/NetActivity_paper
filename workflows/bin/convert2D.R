#! /usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Add distance and CpGs distribution to convert to 2D
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
setPrefix <- args[1]

## Load libraries
library(DelayedMatrixStats)
library(HDF5Array)
library(SummarizedExperiment)
library(rhdf5)

SE <- loadHDF5SummarizedExperiment(dir = "./", prefix = setPrefix)

## Prepare quantile
h5createFile("assay_conv_2D.h5")
h5createDataset("assay_conv_2D.h5", "methy", 
                dims = list(1, nrow(SE), 2, ncol(SE)),
                chunk = c(1, 22736, 2, 1))

dist <-  start(SE)[-1] - start(SE)[-nrow(SE)] 
dist[dist < 0] <- 1e6
dist <- c(dist, 1e6)

for (x in seq_len(ncol(SE))){
  print(x)
  val <- cbind(data.matrix(assay(SE[, x])), dist)
  h5write(val, "assay_conv_2D.h5", "methy", index = list(1, NULL, NULL, x))
} 
h5closeAll()



